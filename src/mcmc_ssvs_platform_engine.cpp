// C++ MCMC engine for the multi-platform SSVS method
// Implements Shaddox et al. (2020, Biostatistics) extension of Peterson et al. (2015).
// Extends the single-platform SSVS (spike-and-slab column-wise Gibbs) with a
// second hierarchical MRF layer that couples group-similarity parameters across
// S data platforms (e.g., metabolomics + transcriptomics).
//
// Phases per iteration:
//   1. Per-platform precision Gibbs (column-wise spike-slab) — identical to SSVS engine
//   2. Per-platform Theta update (MODIFIED: includes platform coupling term via Phi)
//   3. Per-platform nu update — identical to SSVS engine
//   4. Phi update (NEW: platform similarity via MH)
//   5. w_plat update (NEW: platform-level edge log-odds via MH)

#include "wangli_internal.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Helper: stable logistic sigmoid
static inline double stable_sigmoid(double x) {
  if (x > 0.0) {
    double ex = std::exp(-x);
    return 1.0 / (1.0 + ex);
  } else {
    double ex = std::exp(x);
    return ex / (1.0 + ex);
  }
}

// Helper: sum over all edges of [logC(Theta_old, nu_ij) - logC(Theta_new, nu_ij)
//   + 2*(theta_new - theta_old) * g_{k,ij} * g_{m,ij}]
static double sum_logC_diff_edges_plat(const arma::mat& Theta_old, const arma::mat& Theta_new,
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

// Helper: sum over group pairs of [logC(Phi_old, w(k,m)) - logC(Phi_new, w(k,m))
//   + 2*(phi_new - phi_old) * xi[s1](k,m) * xi[s2](k,m)]
// xi_mat is K x K binary: xi_mat(k,m) = 1 if Theta[s](k,m) > 0 for platform s
static double sum_logC_diff_groups(const arma::mat& Phi_old, const arma::mat& Phi_new,
                                    const arma::mat& w_plat,
                                    const std::vector<arma::umat>& xi_mats,
                                    int s1, int s2, double phi_diff, int K) {
  double total = 0.0;
  for (int k = 0; k < K - 1; k++) {
    for (int m = k + 1; m < K; m++) {
      double w_val = w_plat(k, m);
      total += calc_mrf_logC_scalar(Phi_old, w_val)
             - calc_mrf_logC_scalar(Phi_new, w_val)
             + 2.0 * phi_diff * (double)xi_mats[s1](k, m) * (double)xi_mats[s2](k, m);
    }
  }
  return total;
}


// [[Rcpp::export]]
Rcpp::List mcmc_ssvs_platform_engine_cpp(
    Rcpp::List platform_S_lists,
    Rcpp::List platform_n_vecs,
    Rcpp::IntegerVector p_vec_r,
    int K, int burnin, int nsave, int thin,
    double v0, double v1, double lambda_param,
    double a, double b_beta,
    double alpha, double beta_rate,
    double w_prior, double alpha_prop, double beta_prop,
    double a_prop, double b_prop,
    double eta, double kappa,
    double eta_prop, double kappa_prop,
    double d_plat, double f_plat,
    double d_prop, double f_prop,
    double u_prior,
    Rcpp::List Theta_init_list,
    Rcpp::List nu_init_list,
    Rcpp::List C_init_list,
    Rcpp::List Sig_init_list,
    Rcpp::NumericMatrix Phi_init_r,
    Rcpp::NumericMatrix w_init_r,
    bool verbose, int print_every) {

  int S = p_vec_r.size();
  std::vector<int> p_vec(S);
  for (int s = 0; s < S; s++) p_vec[s] = p_vec_r[s];

  // =========================================================================
  // Convert per-platform data from R lists
  // =========================================================================
  std::vector<arma::cube> S_arr(S);
  std::vector<arma::vec> n_vecs(S);
  std::vector<arma::cube> C_arr(S);
  std::vector<arma::cube> Sig_arr(S);
  std::vector<arma::ucube> adj_arr(S);
  std::vector<arma::cube> tau_arr(S);
  std::vector<arma::cube> pii_mat(S);
  std::vector<arma::mat> Theta(S);
  std::vector<arma::mat> nu(S);

  for (int s = 0; s < S; s++) {
    int p_s = p_vec[s];

    // S-matrices
    Rcpp::List S_list_s = Rcpp::as<Rcpp::List>(platform_S_lists[s]);
    S_arr[s].set_size(p_s, p_s, K);
    for (int k = 0; k < K; k++) {
      Rcpp::NumericMatrix Sk = Rcpp::as<Rcpp::NumericMatrix>(S_list_s[k]);
      S_arr[s].slice(k) = arma::mat(Sk.begin(), p_s, p_s, true);
    }

    // Sample sizes
    Rcpp::NumericVector nv = Rcpp::as<Rcpp::NumericVector>(platform_n_vecs[s]);
    n_vecs[s] = arma::vec(nv.begin(), K, true);

    // Theta[s]
    Rcpp::NumericMatrix Th = Rcpp::as<Rcpp::NumericMatrix>(Theta_init_list[s]);
    Theta[s] = arma::mat(Th.begin(), K, K, true);

    // nu[s]
    Rcpp::NumericMatrix Nu = Rcpp::as<Rcpp::NumericMatrix>(nu_init_list[s]);
    nu[s] = arma::mat(Nu.begin(), p_s, p_s, true);

    // C_arr[s] and Sig_arr[s] from flat vectors
    C_arr[s].set_size(p_s, p_s, K);
    Sig_arr[s].set_size(p_s, p_s, K);
    {
      Rcpp::NumericVector c_flat_r = Rcpp::as<Rcpp::NumericVector>(C_init_list[s]);
      Rcpp::NumericVector s_flat_r = Rcpp::as<Rcpp::NumericVector>(Sig_init_list[s]);
      arma::vec c_flat(c_flat_r.begin(), c_flat_r.size(), true);
      arma::vec s_flat(s_flat_r.begin(), s_flat_r.size(), true);
      int idx = 0;
      for (int k = 0; k < K; k++)
        for (int j = 0; j < p_s; j++)
          for (int i = 0; i < p_s; i++) {
            C_arr[s](i, j, k) = c_flat(idx);
            Sig_arr[s](i, j, k) = s_flat(idx);
            idx++;
          }
    }

    // Build initial adjacency from C
    adj_arr[s].set_size(p_s, p_s, K);
    adj_arr[s].zeros();
    for (int k = 0; k < K; k++)
      for (int i = 0; i < p_s; i++)
        for (int j = 0; j < p_s; j++)
          if (i != j && std::abs(C_arr[s](i, j, k)) > 1e-5)
            adj_arr[s](i, j, k) = 1u;

    // tau array
    tau_arr[s].set_size(p_s, p_s, K);
    for (int k = 0; k < K; k++)
      for (int i = 0; i < p_s; i++)
        for (int j = 0; j < p_s; j++)
          tau_arr[s](i, j, k) = (adj_arr[s](i, j, k) == 1u) ? v1 : v0;

    // pii_mat
    pii_mat[s].set_size(p_s, p_s, K);
    for (int k = 0; k < K; k++)
      pii_mat[s].slice(k).fill(2.0 / (double)(p_s - 1));
  }

  // Platform-level initialization
  arma::mat Phi(Phi_init_r.begin(), S, S, true);
  arma::mat w_plat(w_init_r.begin(), K, K, true);

  // xi matrices: xi[s](k,m) = 1 if Theta[s](k,m) > 0
  std::vector<arma::umat> xi_mats(S);
  for (int s = 0; s < S; s++) {
    xi_mats[s].set_size(K, K);
    xi_mats[s].zeros();
    for (int k = 0; k < K; k++)
      for (int m = 0; m < K; m++)
        if (Theta[s](k, m) > 0.0) xi_mats[s](k, m) = 1u;
  }

  // Precompute ind_noi_all per platform
  std::vector<arma::umat> ind_noi_all(S);
  for (int s = 0; s < S; s++) {
    int p_s = p_vec[s];
    ind_noi_all[s].set_size(p_s - 1, p_s);
    for (int i = 0; i < p_s; i++) {
      int idx = 0;
      for (int j = 0; j < p_s; j++) {
        if (j != i) ind_noi_all[s](idx++, i) = (arma::uword)j;
      }
    }
  }

  // =========================================================================
  // Storage
  // =========================================================================
  int nmc = nsave * thin;
  int niter = burnin + nmc;

  // Per-platform saves
  std::vector<std::vector<double>> C_save_vecs(S);
  std::vector<std::vector<double>> Sig_save_vecs(S);
  std::vector<std::vector<int>> adj_save_vecs(S);
  std::vector<arma::cube> Theta_save(S);
  std::vector<arma::cube> nu_save(S);

  for (int s = 0; s < S; s++) {
    int p_s = p_vec[s];
    C_save_vecs[s].assign(p_s * p_s * K * nsave, 0.0);
    Sig_save_vecs[s].assign(p_s * p_s * K * nsave, 0.0);
    adj_save_vecs[s].assign(p_s * p_s * K * nsave, 0);
    Theta_save[s].set_size(K, K, nsave);
    Theta_save[s].zeros();
    nu_save[s].set_size(p_s, p_s, nsave);
    nu_save[s].zeros();
  }

  // Platform-level saves
  arma::cube Phi_save(S, S, nsave, arma::fill::zeros);
  arma::cube w_save(K, K, nsave, arma::fill::zeros);

  // Acceptance tracking — per-platform
  std::vector<arma::mat> ar_gamma(S);
  std::vector<arma::mat> ar_theta(S);
  std::vector<arma::mat> n_within_model(S);
  std::vector<arma::mat> ar_nu(S);
  for (int s = 0; s < S; s++) {
    ar_gamma[s].set_size(K, K); ar_gamma[s].zeros();
    ar_theta[s].set_size(K, K); ar_theta[s].zeros();
    n_within_model[s].set_size(K, K); n_within_model[s].zeros();
    int p_s = p_vec[s];
    ar_nu[s].set_size(p_s, p_s); ar_nu[s].zeros();
  }

  // Platform-level acceptance tracking
  arma::mat ar_phi_between(S, S, arma::fill::zeros);
  arma::mat ar_phi_within(S, S, arma::fill::zeros);
  arma::mat n_phi_within(S, S, arma::fill::zeros);
  arma::mat ar_w(K, K, arma::fill::zeros);

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
    // Phase 1: Update (C_k, Sig_k, adj_k) per platform — column-wise Gibbs
    // -----------------------------------------------------------------------
    for (int s = 0; s < S; s++) {
      int p_s = p_vec[s];
      int pm1 = p_s - 1;

      for (int cur_graph = 0; cur_graph < K; cur_graph++) {
        for (int col_i = 0; col_i < p_s; col_i++) {
          arma::uvec noi = ind_noi_all[s].col(col_i);

          // Extract tau_temp
          arma::vec tau_temp(pm1);
          for (int r = 0; r < pm1; r++)
            tau_temp(r) = tau_arr[s](noi(r), col_i, cur_graph);

          // Sig11, Sig12
          arma::mat Sig11(pm1, pm1);
          arma::vec Sig12(pm1);
          for (int r = 0; r < pm1; r++) {
            Sig12(r) = Sig_arr[s](noi(r), col_i, cur_graph);
            for (int c = 0; c < pm1; c++)
              Sig11(r, c) = Sig_arr[s](noi(r), noi(c), cur_graph);
          }
          double Sig_ii = Sig_arr[s](col_i, col_i, cur_graph);

          // invC11 = Sig11 - Sig12 * Sig12' / Sig_ii
          arma::mat invC11 = Sig11 - Sig12 * Sig12.t() / Sig_ii;
          invC11 = 0.5 * (invC11 + invC11.t());

          // S_noi_i, S_ii
          arma::vec S_noi_i(pm1);
          for (int r = 0; r < pm1; r++)
            S_noi_i(r) = S_arr[s](noi(r), col_i, cur_graph);
          double S_ii = S_arr[s](col_i, col_i, cur_graph);

          // Ci = (S_ii + lambda) * invC11 + diag(1/tau_temp)
          arma::mat Ci = (S_ii + lambda_param) * invC11;
          for (int r = 0; r < pm1; r++)
            Ci(r, r) += 1.0 / tau_temp(r);
          Ci = 0.5 * (Ci + Ci.t());

          // Cholesky
          arma::mat Ci_chol;
          bool chol_ok = safe_chol_try(Ci_chol, Ci, "upper");
          if (!chol_ok) continue;

          // mu_i = -Ci^{-1} * S_noi_i
          arma::vec y = arma::solve(arma::trimatl(Ci_chol.t()), S_noi_i,
                                     arma::solve_opts::fast);
          arma::vec mu_i = -arma::solve(arma::trimatu(Ci_chol), y,
                                         arma::solve_opts::fast);

          // Sample epsilon
          arma::vec z_rand = arma::randn(pm1);
          arma::vec epsilon = mu_i + arma::solve(arma::trimatu(Ci_chol), z_rand,
                                                  arma::solve_opts::fast);

          // Update off-diagonal precision
          for (int r = 0; r < pm1; r++) {
            C_arr[s](noi(r), col_i, cur_graph) = epsilon(r);
            C_arr[s](col_i, noi(r), cur_graph) = epsilon(r);
          }

          // Sample diagonal
          double a_gam = 0.5 * n_vecs[s](cur_graph) + 1.0;
          double b_gam = 0.5 * (S_ii + lambda_param);
          double gam = R::rgamma(a_gam, 1.0 / b_gam);

          double c_val = arma::as_scalar(epsilon.t() * invC11 * epsilon);
          C_arr[s](col_i, col_i, cur_graph) = gam + c_val;

          // Rank-1 covariance update
          arma::vec invC11_eps = invC11 * epsilon;
          arma::mat Sig11_new = invC11 + invC11_eps * invC11_eps.t() / gam;
          for (int r = 0; r < pm1; r++)
            for (int c = 0; c < pm1; c++)
              Sig_arr[s](noi(r), noi(c), cur_graph) = Sig11_new(r, c);

          arma::vec Sig12_new = -invC11_eps / gam;
          for (int r = 0; r < pm1; r++) {
            Sig_arr[s](noi(r), col_i, cur_graph) = Sig12_new(r);
            Sig_arr[s](col_i, noi(r), cur_graph) = Sig12_new(r);
          }
          Sig_arr[s](col_i, col_i, cur_graph) = 1.0 / gam;

          // Update edge indicators via spike-slab posterior
          arma::vec mrf_sum(pm1, arma::fill::zeros);
          for (int m = 0; m < K; m++)
            for (int r = 0; r < pm1; r++)
              mrf_sum(r) += Theta[s](cur_graph, m) * pii_mat[s](noi(r), col_i, m);

          arma::vec nu_noi(pm1);
          for (int r = 0; r < pm1; r++)
            nu_noi(r) = nu[s](noi(r), col_i);
          arma::vec log_odds = nu_noi + 2.0 * mrf_sum;

          arma::vec pii_cur(pm1);
          for (int r = 0; r < pm1; r++)
            pii_cur(r) = stable_sigmoid(log_odds(r));

          for (int r = 0; r < pm1; r++)
            pii_mat[s](noi(r), col_i, cur_graph) = pii_cur(r);

          // Spike-slab weights and Bernoulli draw
          for (int r = 0; r < pm1; r++) {
            double eps_r = epsilon(r);
            double eps2 = eps_r * eps_r;

            double w1 = -0.5 * log_v0 - 0.5 * eps2 / v0 + std::log(1.0 - pii_cur(r) + 1e-300);
            double w2 = -0.5 * log_v1 - 0.5 * eps2 / v1 + std::log(pii_cur(r) + 1e-300);

            double w_max = std::max(w1, w2);
            double prob = std::exp(w2 - w_max) / (std::exp(w1 - w_max) + std::exp(w2 - w_max));

            int z = (R::runif(0.0, 1.0) < prob) ? 1 : 0;

            double v_new = z ? v1 : v0;
            tau_arr[s](noi(r), col_i, cur_graph) = v_new;
            tau_arr[s](col_i, noi(r), cur_graph) = v_new;

            adj_arr[s](noi(r), col_i, cur_graph) = (arma::uword)z;
            adj_arr[s](col_i, noi(r), cur_graph) = (arma::uword)z;
          }

        } // end col_i
      } // end cur_graph
    } // end platform s

    // -----------------------------------------------------------------------
    // Phase 2: Update Theta[s](k,m) per platform — with platform coupling
    // -----------------------------------------------------------------------
    for (int s = 0; s < S; s++) {
      int p_s = p_vec[s];

      for (int k = 0; k < K - 1; k++) {
        for (int m = k + 1; m < K; m++) {

          double theta_curr = Theta[s](k, m);
          bool curr_is_zero = (theta_curr == 0.0);
          bool old_xi = !curr_is_zero;  // xi[s](k,m) = 1 if theta > 0

          // Between-model move
          double theta_prop;
          if (curr_is_zero) {
            theta_prop = R::rgamma(alpha_prop, beta_prop); // scale param
          } else {
            theta_prop = 0.0;
          }

          bool new_xi = (theta_prop > 0.0);

          arma::mat Theta_p = Theta[s];
          Theta_p(k, m) = theta_prop;
          Theta_p(m, k) = theta_prop;

          double sum_edges = sum_logC_diff_edges_plat(Theta[s], Theta_p, nu[s], adj_arr[s],
                                                       k, m, theta_prop - theta_curr, p_s);

          // Platform coupling term: only applies when xi changes and S > 1
          double platform_term = 0.0;
          if (S > 1 && old_xi != new_xi) {
            // Compute: sum over other platforms t of Phi(s,t) * xi[t](k,m)
            // plus the w_plat and u_prior terms from the MRF
            double phi_xi_sum = 0.0;
            for (int t = 0; t < S; t++) {
              if (t == s) continue;
              phi_xi_sum += Phi(s, t) * (double)xi_mats[t](k, m);
            }

            if (new_xi) {
              // 0 -> 1: adding xi[s](k,m)
              platform_term = w_plat(k, m) + 2.0 * phi_xi_sum
                            + std::log(u_prior + 1e-300) - std::log(1.0 - u_prior + 1e-300);
            } else {
              // 1 -> 0: removing xi[s](k,m)
              platform_term = -(w_plat(k, m) + 2.0 * phi_xi_sum)
                            + std::log(1.0 - u_prior + 1e-300) - std::log(u_prior + 1e-300);
            }
          }

          // Log MH ratio
          double log_ar;
          if (theta_prop == 0.0) {
            // Proposing spike from slab
            log_ar = alpha_prop * std::log(beta_prop) - ::lgamma(alpha_prop)
                   + ::lgamma(alpha) - alpha * std::log(beta_rate)
                   - (alpha - alpha_prop) * std::log(theta_curr)
                   + (beta_rate - beta_prop) * theta_curr
                   + sum_edges
                   + std::log(1.0 - w_prior) - std::log(w_prior)
                   + platform_term;
          } else {
            // Proposing slab from spike
            log_ar = alpha * std::log(beta_rate) - ::lgamma(alpha)
                   + ::lgamma(alpha_prop) - alpha_prop * std::log(beta_prop)
                   - (alpha - alpha_prop) * std::log(theta_prop)
                   - (beta_rate - beta_prop) * theta_prop
                   + sum_edges
                   + std::log(w_prior) - std::log(1.0 - w_prior)
                   + platform_term;
          }

          if (std::isfinite(log_ar) && log_ar > std::log(R::runif(0.0, 1.0))) {
            Theta[s](k, m) = theta_prop;
            Theta[s](m, k) = theta_prop;
            // Update xi
            xi_mats[s](k, m) = new_xi ? 1u : 0u;
            xi_mats[s](m, k) = new_xi ? 1u : 0u;
            ar_gamma[s](k, m) += 1.0 / (double)niter;
          }

          // Within-model (slab) move
          if (Theta[s](k, m) != 0.0) {
            double theta_curr2 = Theta[s](k, m);
            double theta_prop2 = R::rgamma(alpha_prop, beta_prop);
            n_within_model[s](k, m) += 1.0;

            arma::mat Theta_p2 = Theta[s];
            Theta_p2(k, m) = theta_prop2;
            Theta_p2(m, k) = theta_prop2;

            double sum_edges2 = sum_logC_diff_edges_plat(Theta[s], Theta_p2, nu[s], adj_arr[s],
                                                          k, m, theta_prop2 - theta_curr2, p_s);

            double log_theta_ar = (alpha - alpha_prop) * (std::log(theta_prop2) - std::log(theta_curr2))
                                + (beta_rate - beta_prop) * (theta_curr2 - theta_prop2)
                                + sum_edges2;

            if (std::isfinite(log_theta_ar) && log_theta_ar > std::log(R::runif(0.0, 1.0))) {
              Theta[s](k, m) = theta_prop2;
              Theta[s](m, k) = theta_prop2;
              ar_theta[s](k, m) += 1.0;
            }
          }

        } // end m
      } // end k
    } // end platform s

    // -----------------------------------------------------------------------
    // Phase 3: Update nu[s] per platform
    // -----------------------------------------------------------------------
    for (int s = 0; s < S; s++) {
      int p_s = p_vec[s];
      for (int ii = 0; ii < p_s - 1; ii++) {
        for (int jj = ii + 1; jj < p_s; jj++) {
          double q_star = R::rbeta(a_prop, b_prop);
          double nu_prop = std::log(q_star) - std::log(1.0 - q_star);

          double sum_g = 0.0;
          for (int k = 0; k < K; k++) sum_g += (double)adj_arr[s](ii, jj, k);

          double log_nu_ar = (nu_prop - nu[s](ii, jj)) * (sum_g + a - a_prop)
            - (a + b_beta - a_prop - b_prop) * std::log(1.0 + std::exp(nu_prop))
            - calc_mrf_logC_scalar(Theta[s], nu_prop)
            + (a + b_beta - a_prop - b_prop) * std::log(1.0 + std::exp(nu[s](ii, jj)))
            + calc_mrf_logC_scalar(Theta[s], nu[s](ii, jj));

          if (log_nu_ar > std::log(R::runif(0.0, 1.0))) {
            nu[s](ii, jj) = nu_prop;
            nu[s](jj, ii) = nu_prop;
            ar_nu[s](ii, jj) += 1.0 / (double)niter;
          }
        }
      }
    } // end platform s

    // -----------------------------------------------------------------------
    // Phase 4: Update Phi(s1, s2) — platform similarity
    // -----------------------------------------------------------------------
    if (S > 1) {
      for (int s1 = 0; s1 < S - 1; s1++) {
        for (int s2 = s1 + 1; s2 < S; s2++) {

          double phi_curr = Phi(s1, s2);

          // Between-model move
          double phi_prop;
          if (phi_curr == 0.0) {
            phi_prop = R::rgamma(eta_prop, kappa_prop); // scale param
          } else {
            phi_prop = 0.0;
          }

          arma::mat Phi_p = Phi;
          Phi_p(s1, s2) = phi_prop;
          Phi_p(s2, s1) = phi_prop;

          // Sum over group pairs
          double sum_groups = sum_logC_diff_groups(Phi, Phi_p, w_plat, xi_mats,
                                                    s1, s2, phi_prop - phi_curr, K);

          double log_ar_phi;
          if (phi_prop == 0.0) {
            // Proposing spike from slab
            log_ar_phi = eta_prop * std::log(kappa_prop) - ::lgamma(eta_prop)
                       + ::lgamma(eta) - eta * std::log(kappa)
                       - (eta - eta_prop) * std::log(phi_curr)
                       + (1.0/kappa - 1.0/kappa_prop) * phi_curr
                       + sum_groups
                       + std::log(1.0 - u_prior) - std::log(u_prior);
          } else {
            // Proposing slab from spike
            log_ar_phi = eta * std::log(1.0/kappa) - ::lgamma(eta)
                       + ::lgamma(eta_prop) - eta_prop * std::log(kappa_prop)
                       - (eta - eta_prop) * std::log(phi_prop)
                       - (1.0/kappa - 1.0/kappa_prop) * phi_prop
                       + sum_groups
                       + std::log(u_prior) - std::log(1.0 - u_prior);
          }

          if (std::isfinite(log_ar_phi) && log_ar_phi > std::log(R::runif(0.0, 1.0))) {
            Phi(s1, s2) = phi_prop;
            Phi(s2, s1) = phi_prop;
            ar_phi_between(s1, s2) += 1.0 / (double)niter;
          }

          // Within-model move
          if (Phi(s1, s2) != 0.0) {
            double phi_curr2 = Phi(s1, s2);
            double phi_prop2 = R::rgamma(eta_prop, kappa_prop);
            n_phi_within(s1, s2) += 1.0;

            arma::mat Phi_p2 = Phi;
            Phi_p2(s1, s2) = phi_prop2;
            Phi_p2(s2, s1) = phi_prop2;

            double sum_groups2 = sum_logC_diff_groups(Phi, Phi_p2, w_plat, xi_mats,
                                                       s1, s2, phi_prop2 - phi_curr2, K);

            double log_phi_ar2 = (eta - eta_prop) * (std::log(phi_prop2) - std::log(phi_curr2))
                               + (1.0/kappa - 1.0/kappa_prop) * (phi_curr2 - phi_prop2)
                               + sum_groups2;

            if (std::isfinite(log_phi_ar2) && log_phi_ar2 > std::log(R::runif(0.0, 1.0))) {
              Phi(s1, s2) = phi_prop2;
              Phi(s2, s1) = phi_prop2;
              ar_phi_within(s1, s2) += 1.0;
            }
          }

        } // end s2
      } // end s1
    } // end if S > 1

    // -----------------------------------------------------------------------
    // Phase 5: Update w_plat(k,m) — platform-level edge log-odds
    // -----------------------------------------------------------------------
    if (S > 1) {
      for (int k = 0; k < K - 1; k++) {
        for (int m = k + 1; m < K; m++) {
          double q_star = R::rbeta(d_prop, f_prop);
          double w_prop = std::log(q_star) - std::log(1.0 - q_star);

          // Sum xi[s](k,m) across platforms
          double sum_xi = 0.0;
          for (int s = 0; s < S; s++) sum_xi += (double)xi_mats[s](k, m);

          double log_w_ar = (w_prop - w_plat(k, m)) * (sum_xi + d_plat - d_prop)
            - (d_plat + f_plat - d_prop - f_prop) * std::log(1.0 + std::exp(w_prop))
            - calc_mrf_logC_scalar(Phi, w_prop)
            + (d_plat + f_plat - d_prop - f_prop) * std::log(1.0 + std::exp(w_plat(k, m)))
            + calc_mrf_logC_scalar(Phi, w_plat(k, m));

          if (log_w_ar > std::log(R::runif(0.0, 1.0))) {
            w_plat(k, m) = w_prop;
            w_plat(m, k) = w_prop;
            ar_w(k, m) += 1.0 / (double)niter;
          }
        }
      }
    }

    // -----------------------------------------------------------------------
    // Store samples
    // -----------------------------------------------------------------------
    if (iter > burnin && ((iter - burnin) % thin == 0)) {
      if (save_idx < nsave) {
        for (int s = 0; s < S; s++) {
          int p_s = p_vec[s];
          int base = save_idx * p_s * p_s * K;
          for (int k = 0; k < K; k++) {
            int off = base + k * p_s * p_s;
            for (int j = 0; j < p_s; j++) {
              for (int i = 0; i < p_s; i++) {
                int idx = off + j * p_s + i;
                C_save_vecs[s][idx] = C_arr[s](i, j, k);
                Sig_save_vecs[s][idx] = Sig_arr[s](i, j, k);
                adj_save_vecs[s][idx] = (int)adj_arr[s](i, j, k);
              }
            }
          }
          Theta_save[s].slice(save_idx) = Theta[s];
          nu_save[s].slice(save_idx) = nu[s];
        }
        Phi_save.slice(save_idx) = Phi;
        w_save.slice(save_idx) = w_plat;
        save_idx++;
      }
    }
  } // end iter

  // =========================================================================
  // Post-process acceptance rates
  // =========================================================================
  for (int s = 0; s < S; s++) {
    for (int k = 0; k < K - 1; k++)
      for (int m = k + 1; m < K; m++)
        if (n_within_model[s](k, m) > 0.0)
          ar_theta[s](k, m) /= n_within_model[s](k, m);
  }
  for (int s1 = 0; s1 < S - 1; s1++)
    for (int s2 = s1 + 1; s2 < S; s2++)
      if (n_phi_within(s1, s2) > 0.0)
        ar_phi_within(s1, s2) /= n_phi_within(s1, s2);

  // =========================================================================
  // Return results
  // =========================================================================
  Rcpp::List platform_results(S);
  for (int s = 0; s < S; s++) {
    int p_s = p_vec[s];
    Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(p_s, p_s, K, nsave);

    Rcpp::NumericVector C_out(C_save_vecs[s].begin(), C_save_vecs[s].end());
    C_out.attr("dim") = dims;

    Rcpp::NumericVector Sig_out(Sig_save_vecs[s].begin(), Sig_save_vecs[s].end());
    Sig_out.attr("dim") = dims;

    Rcpp::IntegerVector adj_out(adj_save_vecs[s].begin(), adj_save_vecs[s].end());
    adj_out.attr("dim") = dims;

    platform_results[s] = Rcpp::List::create(
      Rcpp::_["C_save"] = C_out,
      Rcpp::_["Sig_save"] = Sig_out,
      Rcpp::_["adj_save"] = adj_out,
      Rcpp::_["Theta_save"] = Theta_save[s],
      Rcpp::_["nu_save"] = nu_save[s],
      Rcpp::_["ar_gamma"] = ar_gamma[s],
      Rcpp::_["ar_theta"] = ar_theta[s],
      Rcpp::_["ar_nu"] = ar_nu[s]
    );
  }

  return Rcpp::List::create(
    Rcpp::_["platform_results"] = platform_results,
    Rcpp::_["Phi_save"] = Phi_save,
    Rcpp::_["w_save"] = w_save,
    Rcpp::_["ar_phi_between"] = ar_phi_between,
    Rcpp::_["ar_phi_within"] = ar_phi_within,
    Rcpp::_["ar_w"] = ar_w
  );
}
