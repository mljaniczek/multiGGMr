// C++ MCMC engine for multiGGMr
// Moves the entire per-iteration loop into C++ to avoid R<->C++ overhead.
// All internal functions come from wangli_internal.h.

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
Rcpp::List mcmc_engine_cpp(Rcpp::List S_list_r,
                            Rcpp::NumericVector n_vec_r,
                            int burnin, int nsave, int thin,
                            double b_prior, Rcpp::NumericMatrix D_prior_r,
                            double a, double b_beta,
                            double alpha, double beta_rate,
                            double w_prior,
                            double alpha_prop, double beta_prop,
                            double a_prop, double b_prop,
                            Rcpp::NumericMatrix Theta_init_r,
                            Rcpp::NumericMatrix nu_init_r,
                            Rcpp::NumericVector C_init_r,
                            bool verbose, int print_every) {

  int K = S_list_r.size();
  arma::mat D_prior(D_prior_r.begin(), D_prior_r.nrow(), D_prior_r.ncol(), true);
  int p = (int)D_prior.n_rows;

  // Convert S_list
  arma::cube S_arr(p, p, K);
  for (int k = 0; k < K; k++) {
    Rcpp::NumericMatrix Sk = Rcpp::as<Rcpp::NumericMatrix>(S_list_r[k]);
    S_arr.slice(k) = arma::mat(Sk.begin(), p, p, true);
  }

  arma::vec n_vec(n_vec_r.begin(), K, true);

  // Posterior hyperparameters per group
  arma::vec b_post(K);
  arma::cube D_post(p, p, K);
  for (int k = 0; k < K; k++) {
    b_post(k) = b_prior + n_vec(k);
    D_post.slice(k) = D_prior + S_arr.slice(k);
  }

  // Initialize Theta, nu, C_arr, adj_arr
  arma::mat Theta(Theta_init_r.begin(), K, K, true);
  arma::mat nu(nu_init_r.begin(), p, p, true);

  // C_init comes as a flat vector (p*p*K)
  arma::cube C_arr(p, p, K);
  {
    arma::vec c_flat(C_init_r.begin(), C_init_r.size(), true);
    int idx = 0;
    for (int k = 0; k < K; k++)
      for (int j = 0; j < p; j++)
        for (int i = 0; i < p; i++)
          C_arr(i, j, k) = c_flat(idx++);
  }

  // Build initial adjacency from C
  arma::ucube adj_arr(p, p, K, arma::fill::zeros);
  for (int k = 0; k < K; k++) {
    for (int i = 0; i < p; i++)
      for (int j = 0; j < p; j++)
        if (i != j && std::abs(C_arr(i, j, k)) > 1e-5)
          adj_arr(i, j, k) = 1u;
  }

  // Storage
  int nmc = nsave * thin;
  int niter = burnin + nmc;

  arma::cube C_save_flat(p * p * K, 1, nsave, arma::fill::zeros);
  arma::cube Sig_save_flat(p * p * K, 1, nsave, arma::fill::zeros);
  arma::field<arma::ucube> adj_save_field(nsave);
  arma::cube Theta_save(K, K, nsave, arma::fill::zeros);
  arma::cube nu_save(p, p, nsave, arma::fill::zeros);

  // Use 4D-like storage via flat indexing
  // C_save[p,p,K,nsave], Sig_save[p,p,K,nsave], adj_save[p,p,K,nsave]
  std::vector<double> C_save_vec(p * p * K * nsave, 0.0);
  std::vector<double> Sig_save_vec(p * p * K * nsave, 0.0);
  std::vector<int> adj_save_vec(p * p * K * nsave, 0);

  // Acceptance tracking
  arma::mat ar_gamma(K, K, arma::fill::zeros);
  arma::mat ar_theta(K, K, arma::fill::zeros);
  arma::mat n_within_model(K, K, arma::fill::zeros);
  arma::mat ar_nu(p, p, arma::fill::zeros);

  int save_idx = 0;

  // =========================================================================
  // MCMC main loop
  // =========================================================================
  for (int iter = 1; iter <= niter; iter++) {
    if (verbose && (iter % print_every == 0)) {
      Rcpp::Rcout << "iter = " << iter << " / " << niter << std::endl;
    }
    Rcpp::checkUserInterrupt();

    // -----------------------------------------------------------------------
    // Step 1: Update (G_k, Omega_k) for each group k
    // -----------------------------------------------------------------------
    for (int cur_graph = 0; cur_graph < K; cur_graph++) {
      arma::mat Ck = C_arr.slice(cur_graph);
      arma::umat adjk = adj_arr.slice(cur_graph);

      // Sanitize adjacency
      adjk = (adjk + adjk.t());
      adjk.transform([](arma::uword x) { return (x > 0u) ? 1u : 0u; });
      adjk.diag().zeros();

      for (int ii = 0; ii < p - 1; ii++) {
        for (int jj = ii + 1; jj < p; jj++) {

          // Step 1(a): MRF conditional
          double logH = compute_log_H(b_prior, D_prior, n_vec(cur_graph),
                                       S_arr.slice(cur_graph), Ck, ii, jj);

          double sim_term = nu(ii, jj);
          for (int m = 0; m < K; m++) {
            sim_term += 2.0 * Theta(cur_graph, m) * (double)adj_arr(ii, jj, m);
          }

          double w_logodds = logH - sim_term;

          if (!std::isfinite(w_logodds)) continue;

          double w_prob = stable_sigmoid(w_logodds);
          if (!std::isfinite(w_prob)) continue;

          int current_ij = (int)adjk(ii, jj);
          int propose_ij = (R::runif(0.0, 1.0) < w_prob) ? 1 : 0;

          if (propose_ij != current_ij) {
            // Step 1(b): Wang-Li exchange
            arma::umat adjk_int = adjk;
            // Ensure it's proper integer-like umat
            GWishartResult out_prior = gwishart_NOij_gibbs_internal(
              b_prior, D_prior, adjk_int, Ck,
              ii, jj, propose_ij, 0, 1);

            if (!out_prior.ok) continue;
            arma::mat C_prop = out_prior.C;

            double r2 = compute_log_GWishart_NOij_pdf(b_prior, D_prior,
                                                       C_prop, ii, jj, current_ij) -
                         compute_log_GWishart_NOij_pdf(b_prior, D_prior,
                                                       C_prop, ii, jj, propose_ij);

            if (std::log(R::runif(0.0, 1.0)) < r2) {
              adjk(ii, jj) = (arma::uword)propose_ij;
              adjk(jj, ii) = (arma::uword)propose_ij;
              current_ij = propose_ij;
            }
          }

          // Step 2(c): Posterior Gibbs update
          arma::umat adjk_int2 = adjk;
          GWishartResult out_g = gwishart_NOij_gibbs_internal(
            b_post(cur_graph), D_post.slice(cur_graph),
            adjk_int2, Ck, ii, jj, current_ij, 0, 0);

          if (!out_g.ok) continue;
          Ck = out_g.C;
          adjk = out_g.adj;
        }
      }

      // Step 3: BIPS clique update
      GWishartResult out_bips = gwishart_BIPS_internal(
        b_post(cur_graph), D_post.slice(cur_graph),
        adjk, Ck, 0, 1);

      if (out_bips.ok) {
        C_arr.slice(cur_graph) = out_bips.C;
      }
      adj_arr.slice(cur_graph) = adjk;
    }

    // -----------------------------------------------------------------------
    // Step 2: Update (theta_km, gamma_km)
    // -----------------------------------------------------------------------
    for (int k = 0; k < K - 1; k++) {
      for (int m = k + 1; m < K; m++) {

        double theta_curr = Theta(k, m);

        // Between-model move
        double theta_prop;
        if (theta_curr == 0.0) {
          theta_prop = R::rgamma(alpha_prop, 1.0 / beta_prop);
        } else {
          theta_prop = 0.0;
        }

        arma::mat Theta_p = Theta;
        Theta_p(k, m) = theta_prop;
        Theta_p(m, k) = theta_prop;

        double sum_edges = sum_logC_diff_edges(Theta, Theta_p, nu, adj_arr,
                                                k, m, theta_prop - theta_curr, p);

        // Log prior ratio
        double log_prior_prop, log_prior_curr;
        if (theta_prop == 0.0) {
          log_prior_prop = std::log(1.0 - w_prior);
        } else {
          log_prior_prop = std::log(w_prior) + R::dgamma(theta_prop, alpha, 1.0 / beta_rate, 1);
        }
        if (theta_curr == 0.0) {
          log_prior_curr = std::log(1.0 - w_prior);
        } else {
          log_prior_curr = std::log(w_prior) + R::dgamma(theta_curr, alpha, 1.0 / beta_rate, 1);
        }

        // Proposal ratio
        double log_q_prop_given_curr, log_q_curr_given_prop;
        if (theta_curr == 0.0) {
          log_q_prop_given_curr = R::dgamma(theta_prop, alpha_prop, 1.0 / beta_prop, 1);
        } else {
          log_q_prop_given_curr = 0.0;
        }
        if (theta_prop == 0.0) {
          log_q_curr_given_prop = R::dgamma(theta_curr, alpha_prop, 1.0 / beta_prop, 1);
        } else {
          log_q_curr_given_prop = 0.0;
        }

        double log_ar = (log_prior_prop - log_prior_curr) + sum_edges +
                         (log_q_curr_given_prop - log_q_prop_given_curr);

        if (std::isfinite(log_ar) && log_ar > std::log(R::runif(0.0, 1.0))) {
          Theta(k, m) = theta_prop;
          Theta(m, k) = theta_prop;
          ar_gamma(k, m) += 1.0 / (double)niter;
        }

        // Within-model (slab) move
        if (Theta(k, m) != 0.0) {
          double theta_curr2 = Theta(k, m);
          double theta_prop2 = R::rgamma(alpha_prop, 1.0 / beta_prop);
          n_within_model(k, m) += 1.0;

          arma::mat Theta_p2 = Theta;
          Theta_p2(k, m) = theta_prop2;
          Theta_p2(m, k) = theta_prop2;

          double sum_edges2 = sum_logC_diff_edges(Theta, Theta_p2, nu, adj_arr,
                                                   k, m, theta_prop2 - theta_curr2, p);

          double log_target = (R::dgamma(theta_prop2, alpha, 1.0 / beta_rate, 1) -
                                R::dgamma(theta_curr2, alpha, 1.0 / beta_rate, 1)) +
                               sum_edges2;

          double log_hastings = R::dgamma(theta_curr2, alpha_prop, 1.0 / beta_prop, 1) -
                                 R::dgamma(theta_prop2, alpha_prop, 1.0 / beta_prop, 1);

          double log_theta_ar = log_target + log_hastings;

          if (std::isfinite(log_theta_ar) && log_theta_ar > std::log(R::runif(0.0, 1.0))) {
            Theta(k, m) = theta_prop2;
            Theta(m, k) = theta_prop2;
            ar_theta(k, m) += 1.0;
          }
        }
      }
    }

    // -----------------------------------------------------------------------
    // Step 3: Update nu_ij
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
          arma::mat Sig_k;
          bool inv_ok = safe_inv_sympd_try(Sig_k, C_arr.slice(k));

          int off = base + k * p * p;
          for (int j = 0; j < p; j++) {
            for (int i = 0; i < p; i++) {
              int idx = off + j * p + i;
              C_save_vec[idx] = C_arr(i, j, k);
              Sig_save_vec[idx] = inv_ok ? Sig_k(i, j) : 0.0;
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
  // Wrap flat vectors into R arrays with proper dimensions
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
