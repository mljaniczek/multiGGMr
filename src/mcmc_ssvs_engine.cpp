// C++ MCMC engine for the SSVS (spike-and-slab) method
// Implements Shaddox et al. (2018) scalable multiGGM approach.
// Replaces the G-Wishart + Wang-Li exchange algorithm with column-wise
// Gibbs sampling using spike-and-slab normal priors on precision elements.
// The MRF coupling, Theta, and nu updates are identical to the G-Wishart engine.

#include "wangli_internal.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Helper: stable logistic sigmoid
static inline double stable_sigmoid(double x) {
  if (x > 0.0) {
    double ex = std::exp(-x);
    return ex / (1.0 + ex);
  } else {
    double ex = std::exp(x);
    return 1.0 / (1.0 + ex);
  }
}

// Helper: sum over all edges of [logC(Theta_old, nu_ij) - logC(Theta_new, nu_ij)
//   + 2*(theta_new - theta_old) * g_{k,ij} * g_{m,ij}]
static double sum_logC_diff_edges(const arma::mat& Theta_old, const arma::mat& Theta_new,
                                   const arma::mat& nu, const arma::ucube& adj_arr,
                                   int k, int m, double theta_diff, int p) {
  double total = 0.0;
  for (int ii = 0; ii < p - 1; ii++) {
    for (int jj = ii + 1; jj < p; jj++) {
      double nu_val = nu(ii, jj);
      total += calc_mrf_logC_scalar(Theta_old, nu_val)
             - calc_mrf_logC_scalar(Theta_new, nu_val)
             + 2.0 * theta_diff * (double)adj_arr(ii, jj, k) * (double)adj_arr(ii, jj, m);
    }
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::List mcmc_ssvs_engine_cpp(Rcpp::List S_list_r,
                                 Rcpp::NumericVector n_vec_r,
                                 int burnin, int nsave, int thin,
                                 double v0, double v1, double lambda_param,
                                 double a, double b_beta,
                                 double alpha, double beta_rate,
                                 double w_prior,
                                 double alpha_prop, double beta_prop,
                                 double a_prop, double b_prop,
                                 Rcpp::NumericMatrix Theta_init_r,
                                 Rcpp::NumericMatrix nu_init_r,
                                 Rcpp::NumericVector C_init_r,
                                 Rcpp::NumericVector Sig_init_r,
                                 bool verbose, int print_every) {

  int K = S_list_r.size();
  int p = Rcpp::as<Rcpp::NumericMatrix>(S_list_r[0]).nrow();

  // Convert S_list
  arma::cube S_arr(p, p, K);
  for (int k = 0; k < K; k++) {
    Rcpp::NumericMatrix Sk = Rcpp::as<Rcpp::NumericMatrix>(S_list_r[k]);
    S_arr.slice(k) = arma::mat(Sk.begin(), p, p, true);
  }

  arma::vec n_vec(n_vec_r.begin(), K, true);

  // Initialize Theta, nu
  arma::mat Theta(Theta_init_r.begin(), K, K, true);
  arma::mat nu(nu_init_r.begin(), p, p, true);

  // C_init and Sig_init come as flat vectors (p*p*K)
  arma::cube C_arr(p, p, K);
  arma::cube Sig_arr(p, p, K);
  {
    arma::vec c_flat(C_init_r.begin(), C_init_r.size(), true);
    arma::vec s_flat(Sig_init_r.begin(), Sig_init_r.size(), true);
    int idx = 0;
    for (int k = 0; k < K; k++)
      for (int j = 0; j < p; j++)
        for (int i = 0; i < p; i++) {
          C_arr(i, j, k) = c_flat(idx);
          Sig_arr(i, j, k) = s_flat(idx);
          idx++;
        }
  }

  // Build initial adjacency from C
  arma::ucube adj_arr(p, p, K, arma::fill::zeros);
  for (int k = 0; k < K; k++) {
    for (int i = 0; i < p; i++)
      for (int j = 0; j < p; j++)
        if (i != j && std::abs(C_arr(i, j, k)) > 1e-5)
          adj_arr(i, j, k) = 1u;
  }

  // tau array: spike or slab variance per element per group
  arma::cube tau_arr(p, p, K);
  for (int k = 0; k < K; k++) {
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        tau_arr(i, j, k) = (adj_arr(i, j, k) == 1u) ? v1 : v0;
      }
    }
  }

  // pii_mat: edge inclusion probabilities per group (initialized to prior)
  // Matlab: pii_mat(:,:,i) = pii(i) where pii = 2/(p-1) typically
  arma::cube pii_mat(p, p, K);
  for (int k = 0; k < K; k++) {
    pii_mat.slice(k).fill(2.0 / (double)(p - 1));
  }

  // Precompute ind_noi_all: for column i, indices of all other columns
  // ind_noi_all is (p-1) x p, each column contains {0,...,p-1}\{i}
  arma::umat ind_noi_all(p - 1, p);
  for (int i = 0; i < p; i++) {
    int idx = 0;
    for (int j = 0; j < p; j++) {
      if (j != i) ind_noi_all(idx++, i) = (arma::uword)j;
    }
  }

  // Storage
  int nmc = nsave * thin;
  int niter = burnin + nmc;

  std::vector<double> C_save_vec(p * p * K * nsave, 0.0);
  std::vector<double> Sig_save_vec(p * p * K * nsave, 0.0);
  std::vector<int> adj_save_vec(p * p * K * nsave, 0);
  arma::cube Theta_save(K, K, nsave, arma::fill::zeros);
  arma::cube nu_save(p, p, nsave, arma::fill::zeros);

  // Acceptance tracking
  arma::mat ar_gamma(K, K, arma::fill::zeros);
  arma::mat ar_theta(K, K, arma::fill::zeros);
  arma::mat n_within_model(K, K, arma::fill::zeros);
  arma::mat ar_nu(p, p, arma::fill::zeros);

  int save_idx = 0;

  double log_v0 = std::log(v0);
  double log_v1 = std::log(v1);

  // =========================================================================
  // MCMC main loop
  // =========================================================================
  for (int iter = 1; iter <= niter; iter++) {
    if (verbose && (iter % print_every == 0)) {
      Rcpp::Rcout << "iter = " << iter << " / " << niter << std::endl;
    }
    Rcpp::checkUserInterrupt();

    // -----------------------------------------------------------------------
    // Phase 1: Update (C_k, Sig_k, adj_k) via column-wise Gibbs (SSVS)
    // -----------------------------------------------------------------------
    for (int cur_graph = 0; cur_graph < K; cur_graph++) {

      for (int col_i = 0; col_i < p; col_i++) {
        // Get indices of all rows except col_i
        arma::uvec noi = ind_noi_all.col(col_i);
        int pm1 = p - 1;

        // Extract sub-matrices
        arma::vec tau_temp(pm1);
        for (int r = 0; r < pm1; r++) {
          tau_temp(r) = tau_arr(noi(r), col_i, cur_graph);
        }

        // Sig11 = Sig(noi, noi, k), Sig12 = Sig(noi, col_i, k)
        arma::mat Sig11(pm1, pm1);
        arma::vec Sig12(pm1);
        for (int r = 0; r < pm1; r++) {
          Sig12(r) = Sig_arr(noi(r), col_i, cur_graph);
          for (int c = 0; c < pm1; c++) {
            Sig11(r, c) = Sig_arr(noi(r), noi(c), cur_graph);
          }
        }
        double Sig_ii = Sig_arr(col_i, col_i, cur_graph);

        // invC11 = Sig11 - Sig12 * Sig12' / Sig(i,i,k)
        arma::mat invC11 = Sig11 - Sig12 * Sig12.t() / Sig_ii;
        invC11 = 0.5 * (invC11 + invC11.t()); // symmetrize

        // S(noi, col_i, k) and S(col_i, col_i, k)
        arma::vec S_noi_i(pm1);
        for (int r = 0; r < pm1; r++) {
          S_noi_i(r) = S_arr(noi(r), col_i, cur_graph);
        }
        double S_ii = S_arr(col_i, col_i, cur_graph);

        // Ci = (S(i,i) + lambda) * invC11 + diag(1/tau_temp)
        arma::mat Ci = (S_ii + lambda_param) * invC11;
        for (int r = 0; r < pm1; r++) {
          Ci(r, r) += 1.0 / tau_temp(r);
        }
        Ci = 0.5 * (Ci + Ci.t()); // symmetrize

        // Cholesky decomposition (upper triangular for solving)
        arma::mat Ci_chol;
        bool chol_ok = safe_chol_try(Ci_chol, Ci, "upper");
        if (!chol_ok) continue; // skip this column if Cholesky fails

        // mu_i = -Ci^{-1} * S(noi, i)
        // Using Ci_chol (upper): solve Ci_chol' * y = S_noi_i, then Ci_chol * mu_i = -y
        arma::vec y = arma::solve(arma::trimatl(Ci_chol.t()), S_noi_i,
                                   arma::solve_opts::fast);
        arma::vec mu_i = -arma::solve(arma::trimatu(Ci_chol), y,
                                       arma::solve_opts::fast);

        // Sample epsilon = mu_i + Ci_chol^{-1} * randn(p-1)
        arma::vec z_rand = arma::randn(pm1);
        arma::vec epsilon = mu_i + arma::solve(arma::trimatu(Ci_chol), z_rand,
                                                arma::solve_opts::fast);

        // Update off-diagonal of precision matrix column
        for (int r = 0; r < pm1; r++) {
          C_arr(noi(r), col_i, cur_graph) = epsilon(r);
          C_arr(col_i, noi(r), cur_graph) = epsilon(r);
        }

        // Sample diagonal: gam ~ Gamma(n_k/2 + 1, scale = 2/(S_ii + lambda))
        double a_gam = 0.5 * n_vec(cur_graph) + 1.0;
        double b_gam = 0.5 * (S_ii + lambda_param);
        double gam = R::rgamma(a_gam, 1.0 / b_gam); // R::rgamma uses scale

        // c = epsilon' * invC11 * epsilon
        double c_val = arma::as_scalar(epsilon.t() * invC11 * epsilon);
        C_arr(col_i, col_i, cur_graph) = gam + c_val;

        // Rank-1 covariance update
        arma::vec invC11_eps = invC11 * epsilon;

        // Sig(noi, noi, k) = invC11 + invC11_eps * invC11_eps' / gam
        arma::mat Sig11_new = invC11 + invC11_eps * invC11_eps.t() / gam;
        for (int r = 0; r < pm1; r++) {
          for (int c = 0; c < pm1; c++) {
            Sig_arr(noi(r), noi(c), cur_graph) = Sig11_new(r, c);
          }
        }

        // Sig(noi, i) = -invC11_eps / gam
        arma::vec Sig12_new = -invC11_eps / gam;
        for (int r = 0; r < pm1; r++) {
          Sig_arr(noi(r), col_i, cur_graph) = Sig12_new(r);
          Sig_arr(col_i, noi(r), cur_graph) = Sig12_new(r);
        }
        Sig_arr(col_i, col_i, cur_graph) = 1.0 / gam;

        // Update edge indicators via spike-slab posterior
        // Compute MRF sum: sum_m Theta(cur_graph, m) * pii_mat(noi, col_i, m)
        arma::vec mrf_sum(pm1, arma::fill::zeros);
        for (int m = 0; m < K; m++) {
          for (int r = 0; r < pm1; r++) {
            mrf_sum(r) += Theta(cur_graph, m) * pii_mat(noi(r), col_i, m);
          }
        }

        // pii_mat(noi, col_i, cur_graph) = sigmoid(nu(noi, col_i) + 2*mrf_sum)
        arma::vec nu_noi(pm1);
        for (int r = 0; r < pm1; r++) {
          nu_noi(r) = nu(noi(r), col_i);
        }
        arma::vec log_odds = nu_noi + 2.0 * mrf_sum;

        arma::vec pii_cur(pm1);
        for (int r = 0; r < pm1; r++) {
          // sigmoid: 1 / (1 + exp(-x))
          double x = log_odds(r);
          if (x > 0.0) {
            double ex = std::exp(-x);
            pii_cur(r) = 1.0 / (1.0 + ex);
          } else {
            double ex = std::exp(x);
            pii_cur(r) = ex / (1.0 + ex);
          }
        }

        // Store updated pii_mat for this group
        for (int r = 0; r < pm1; r++) {
          pii_mat(noi(r), col_i, cur_graph) = pii_cur(r);
        }

        // Spike-slab weights
        // w1 = -0.5*log(v0) - 0.5*epsilon^2/v0 + log(1 - pii)
        // w2 = -0.5*log(v1) - 0.5*epsilon^2/v1 + log(pii)
        for (int r = 0; r < pm1; r++) {
          double eps_r = epsilon(r);
          double eps2 = eps_r * eps_r;

          double w1 = -0.5 * log_v0 - 0.5 * eps2 / v0 + std::log(1.0 - pii_cur(r) + 1e-300);
          double w2 = -0.5 * log_v1 - 0.5 * eps2 / v1 + std::log(pii_cur(r) + 1e-300);

          // Numerically stable softmax
          double w_max = std::max(w1, w2);
          double prob = std::exp(w2 - w_max) / (std::exp(w1 - w_max) + std::exp(w2 - w_max));

          // Bernoulli draw
          int z = (R::runif(0.0, 1.0) < prob) ? 1 : 0;

          double v_new = z ? v1 : v0;
          tau_arr(noi(r), col_i, cur_graph) = v_new;
          tau_arr(col_i, noi(r), cur_graph) = v_new;

          adj_arr(noi(r), col_i, cur_graph) = (arma::uword)z;
          adj_arr(col_i, noi(r), cur_graph) = (arma::uword)z;
        }

      } // end col_i loop
    } // end cur_graph loop

    // -----------------------------------------------------------------------
    // Phase 2: Update (theta_km, gamma_km) — identical to G-Wishart engine
    // -----------------------------------------------------------------------
    for (int k = 0; k < K - 1; k++) {
      for (int m = k + 1; m < K; m++) {

        double theta_curr = Theta(k, m);

        // Between-model move
        double theta_prop;
        if (theta_curr == 0.0) {
          // Matlab: gamrnd(alpha_prop, beta_prop) uses scale parameterization
          theta_prop = R::rgamma(alpha_prop, beta_prop);
        } else {
          theta_prop = 0.0;
        }

        arma::mat Theta_p = Theta;
        Theta_p(k, m) = theta_prop;
        Theta_p(m, k) = theta_prop;

        double sum_edges = sum_logC_diff_edges(Theta, Theta_p, nu, adj_arr,
                                                k, m, theta_prop - theta_curr, p);

        // Log MH ratio (matching Matlab's formulation)
        double log_ar;
        if (theta_prop == 0.0) {
          // Proposing spike (theta_prop = 0) from slab (theta_curr > 0)
          // Matlab: log_ar = alpha_prop*log(beta_prop) - log(gamma(alpha_prop))
          //   + log(gamma(alpha)) - alpha*log(beta)
          //   - (alpha - alpha_prop)*log(theta_curr)
          //   + (beta - beta_prop)*theta_curr + sum_edges
          //   + log(1-w) - log(w)
          log_ar = alpha_prop * std::log(beta_prop) - ::lgamma(alpha_prop)
                 + ::lgamma(alpha) - alpha * std::log(beta_rate)
                 - (alpha - alpha_prop) * std::log(theta_curr)
                 + (beta_rate - beta_prop) * theta_curr
                 + sum_edges
                 + std::log(1.0 - w_prior) - std::log(w_prior);
        } else {
          // Proposing slab (theta_prop > 0) from spike (theta_curr = 0)
          log_ar = alpha * std::log(beta_rate) - ::lgamma(alpha)
                 + ::lgamma(alpha_prop) - alpha_prop * std::log(beta_prop)
                 - (alpha - alpha_prop) * std::log(theta_prop)
                 - (beta_rate - beta_prop) * theta_prop
                 + sum_edges
                 + std::log(w_prior) - std::log(1.0 - w_prior);
        }

        if (std::isfinite(log_ar) && log_ar > std::log(R::runif(0.0, 1.0))) {
          Theta(k, m) = theta_prop;
          Theta(m, k) = theta_prop;
          ar_gamma(k, m) += 1.0 / (double)niter;
        }

        // Within-model (slab) move
        if (Theta(k, m) != 0.0) {
          double theta_curr2 = Theta(k, m);
          double theta_prop2 = R::rgamma(alpha_prop, beta_prop);
          n_within_model(k, m) += 1.0;

          arma::mat Theta_p2 = Theta;
          Theta_p2(k, m) = theta_prop2;
          Theta_p2(m, k) = theta_prop2;

          double sum_edges2 = sum_logC_diff_edges(Theta, Theta_p2, nu, adj_arr,
                                                   k, m, theta_prop2 - theta_curr2, p);

          // Matlab: (alpha - alpha_prop)*(log(theta_prop2) - log(theta_curr2))
          //       + (beta - beta_prop)*(theta_curr2 - theta_prop2) + sum_edges
          double log_theta_ar = (alpha - alpha_prop) * (std::log(theta_prop2) - std::log(theta_curr2))
                              + (beta_rate - beta_prop) * (theta_curr2 - theta_prop2)
                              + sum_edges2;

          if (std::isfinite(log_theta_ar) && log_theta_ar > std::log(R::runif(0.0, 1.0))) {
            Theta(k, m) = theta_prop2;
            Theta(m, k) = theta_prop2;
            ar_theta(k, m) += 1.0;
          }
        }
      }
    }

    // -----------------------------------------------------------------------
    // Phase 3: Update nu_ij — identical to G-Wishart engine
    // -----------------------------------------------------------------------
    for (int ii = 0; ii < p - 1; ii++) {
      for (int jj = ii + 1; jj < p; jj++) {
        double q_star = R::rbeta(a_prop, b_prop);
        double nu_prop = std::log(q_star) - std::log(1.0 - q_star);

        double sum_g = 0.0;
        for (int k = 0; k < K; k++) sum_g += (double)adj_arr(ii, jj, k);

        double log_nu_ar = (nu_prop - nu(ii, jj)) * (sum_g + a - a_prop)
          - (a + b_beta - a_prop - b_prop) * std::log(1.0 + std::exp(nu_prop))
          - calc_mrf_logC_scalar(Theta, nu_prop)
          + (a + b_beta - a_prop - b_prop) * std::log(1.0 + std::exp(nu(ii, jj)))
          + calc_mrf_logC_scalar(Theta, nu(ii, jj));

        if (log_nu_ar > std::log(R::runif(0.0, 1.0))) {
          nu(ii, jj) = nu_prop;
          nu(jj, ii) = nu_prop;
          ar_nu(ii, jj) += 1.0 / (double)niter;
        }
      }
    }

    // -----------------------------------------------------------------------
    // Store samples
    // -----------------------------------------------------------------------
    if (iter > burnin && ((iter - burnin) % thin == 0)) {
      if (save_idx < nsave) {
        int base = save_idx * p * p * K;
        for (int k = 0; k < K; k++) {
          int off = base + k * p * p;
          for (int j = 0; j < p; j++) {
            for (int i = 0; i < p; i++) {
              int idx = off + j * p + i;
              C_save_vec[idx] = C_arr(i, j, k);
              Sig_save_vec[idx] = Sig_arr(i, j, k);
              adj_save_vec[idx] = (int)adj_arr(i, j, k);
            }
          }
        }
        Theta_save.slice(save_idx) = Theta;
        nu_save.slice(save_idx) = nu;
        save_idx++;
      }
    }
  }

  // Post-process acceptance rates
  for (int k = 0; k < K - 1; k++) {
    for (int m = k + 1; m < K; m++) {
      if (n_within_model(k, m) > 0.0) {
        ar_theta(k, m) /= n_within_model(k, m);
      }
    }
  }

  // Return results
  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(p, p, K, nsave);

  Rcpp::NumericVector C_out(C_save_vec.begin(), C_save_vec.end());
  C_out.attr("dim") = dims;

  Rcpp::NumericVector Sig_out(Sig_save_vec.begin(), Sig_save_vec.end());
  Sig_out.attr("dim") = dims;

  Rcpp::IntegerVector adj_out(adj_save_vec.begin(), adj_save_vec.end());
  adj_out.attr("dim") = dims;

  return Rcpp::List::create(
    Rcpp::_["C_save"] = C_out,
    Rcpp::_["Sig_save"] = Sig_out,
    Rcpp::_["adj_save"] = adj_out,
    Rcpp::_["Theta_save"] = Theta_save,
    Rcpp::_["nu_save"] = nu_save,
    Rcpp::_["ar_gamma"] = ar_gamma,
    Rcpp::_["ar_theta"] = ar_theta,
    Rcpp::_["ar_nu"] = ar_nu
  );
}
