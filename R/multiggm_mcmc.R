#' Fit multi-graph Gaussian graphical model via MCMC
#'
#' Implements Bayesian inference of multiple Gaussian graphical models following
#' Peterson, Stingo & Vannucci (2015, JASA). Supports both raw data and
#' pre-computed covariance matrices, with optional parallel multi-chain execution.
#'
#' @param data_list Optional list of K data matrices (each n_k x p). If provided,
#'   \code{S_list} and \code{n_vec} are computed automatically. Column-centers
#'   data internally. Exactly one of \code{data_list} or \code{S_list} must be given.
#' @param S_list Optional list of K cross-product matrices (each p x p), where
#'   \code{S_list[[k]] = t(X_k) \%*\% X_k} after column-centering.
#' @param n_vec Integer vector of length K with sample sizes. Required when
#'   \code{S_list} is provided; ignored when \code{data_list} is provided.
#' @param burnin Integer burn-in iterations.
#' @param nsave Integer number of saved posterior draws per chain (after thinning).
#' @param thin Integer thinning interval.
#' @param nchains Number of independent MCMC chains.
#' @param parallel Logical; if TRUE, run chains in parallel using \pkg{parallel}.
#' @param ncores Number of cores/workers to use when \code{parallel=TRUE}.
#' @param seed Optional integer seed. If provided, chains use deterministic offsets.
#' @param hyper Optional named list of hyperparameters. See Details.
#' @param ... Additional arguments forwarded to the single-chain engine.
#'
#' @details
#' The \code{hyper} list may contain any of the following (with defaults from
#' Peterson et al. 2015, Section 5.1):
#' \describe{
#'   \item{\code{b_prior}}{G-Wishart degrees of freedom (default 3).}
#'   \item{\code{D_prior}}{G-Wishart scale matrix, p x p (default \code{diag(p)}).}
#'   \item{\code{a, b}}{Beta(a, b) prior on q_ij = logit^{-1}(nu_ij). Default a=1, b=4
#'     giving prior edge probability ~0.20.}
#'   \item{\code{alpha, beta}}{Gamma(shape=alpha, rate=beta) prior on theta_km.
#'     NOTE: \code{beta} is a RATE parameter. Default alpha=2, beta=5 giving mean=0.4.}
#'   \item{\code{w}}{Bernoulli prior probability that gamma_km=1. Default 0.9.}
#'   \item{\code{alpha_prop, beta_prop}}{Gamma(shape=alpha_prop, rate=beta_prop)
#'     proposal for theta_km. Default alpha_prop=2, beta_prop=5.}
#' }
#'
#' @return An object of class \code{"multiggm_fit"} (single chain) or
#'   \code{"multiggm_fit_list"} (multiple chains).
#' @export
multiggm_mcmc <- function(data_list = NULL, S_list = NULL, n_vec = NULL,
                          burnin = 5e3, nsave = 1e3, thin = 1,
                          nchains = 1,
                          parallel = FALSE, ncores = max(1L, parallel::detectCores() - 1L),
                          seed = NULL,
                          hyper = NULL,
                          ...) {

  # --- Input validation: exactly one of data_list or S_list ---
  has_data <- !is.null(data_list)
  has_S    <- !is.null(S_list)
  if (has_data == has_S) stop("Provide exactly one of 'data_list' or 'S_list'.")

  if (has_data) {
    if (!is.list(data_list) || length(data_list) < 1L)
      stop("data_list must be a non-empty list of data matrices.")
    K <- length(data_list)
    p <- ncol(data_list[[1L]])
    n_vec <- integer(K)
    S_list <- vector("list", K)
    for (k in seq_len(K)) {
      Xk <- as.matrix(data_list[[k]])
      if (ncol(Xk) != p) stop("All data matrices must have the same number of columns.")
      n_vec[k] <- nrow(Xk)
      # Column-center and compute cross-product (not cov(); paper uses X'X)
      Xk <- scale(Xk, center = TRUE, scale = FALSE)
      S_list[[k]] <- crossprod(Xk)
    }
  } else {
    if (!is.list(S_list) || length(S_list) < 1L)
      stop("S_list must be a non-empty list of cross-product matrices.")
    K <- length(S_list)
    n_vec <- as.integer(n_vec)
    if (length(n_vec) != K) stop("n_vec must have length equal to length(S_list).")
  }

  if (!is.null(seed)) seed <- as.integer(seed)

  chain_fun <- function(chain_id) {
    if (!is.null(seed)) set.seed(seed + 1000L * (chain_id - 1L))
    .multiggm_mcmc_single(S_list = S_list, n_vec = n_vec,
                         burnin = burnin, nsave = nsave, thin = thin,
                         chain_id = chain_id, hyper = hyper, ...)
  }

  if (nchains == 1L) {
    fit <- chain_fun(1L)
    class(fit) <- c("multiggm_fit", class(fit))
    return(fit)
  }

  if (parallel) {
    ncores <- as.integer(ncores)
    if (.Platform$OS.type == "unix") {
      fits <- parallel::mclapply(seq_len(nchains), chain_fun, mc.cores = ncores)
    } else {
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, varlist = c("S_list", "n_vec", "burnin", "nsave",
                                               "thin", "seed", "hyper"),
                              envir = environment())
      fits <- parallel::parLapply(cl, seq_len(nchains), chain_fun)
    }
  } else {
    fits <- lapply(seq_len(nchains), chain_fun)
  }

  out <- list(chains = fits,
              call = match.call(),
              K = K,
              nchains = nchains)
  class(out) <- "multiggm_fit_list"
  out
}

#' Internal single-chain MCMC engine
#'
#' Implements the MCMC sampler from Peterson et al. (2015) Appendix A using the
#' Wang-Li (2012) G-Wishart exchange algorithm (ported to RcppArmadillo).
#'
#' @keywords internal
.multiggm_mcmc_single <- function(S_list, n_vec, burnin, nsave, thin,
                                 chain_id = 1L,
                                 Theta_init = NULL,
                                 nu_init = NULL,
                                 C_init = NULL,
                                 disp = FALSE,
                                 hyper = NULL, ...) {

  stopifnot(is.list(S_list), length(S_list) >= 1L)
  K <- length(S_list)
  p <- nrow(S_list[[1L]])

  # --- Default hyperparameters (Peterson et al. 2015, Section 5.1) ---
  if (is.null(hyper)) {
    hyper <- list(
      b_prior = 3,                # G-Wishart df (Section 3.5)
      D_prior = diag(p),          # G-Wishart scale (Section 3.5)
      a = 1, b = 4,               # Beta(1,4) prior on q_ij -> prior P(edge) ~ 0.20 (Section 3.4)
      alpha = 2, beta = 5,        # Gamma(2, rate=5) prior on theta_km -> mean=0.4 (Section 3.3)
      w = 0.9,                    # P(gamma_km=1) -> strong relatedness belief (Section 3.3)
      alpha_prop = 2, beta_prop = 5  # Gamma proposal for theta (Appendix A.2)
    )
  }

  # Unpack hyperparameters
  b_prior    <- hyper$b_prior
  D_prior    <- hyper$D_prior
  a          <- hyper$a
  b          <- hyper$b
  alpha      <- hyper$alpha
  beta       <- hyper$beta       # NOTE: this is a RATE parameter
  my_w       <- hyper$w
  alpha_prop <- hyper$alpha_prop
  beta_prop  <- hyper$beta_prop  # NOTE: this is a RATE parameter

  # Beta proposal for nu (Appendix A.3: "propose q* from Beta(2,4)")
  a_prop <- 2
  b_prop <- 4

  # --- Initialization ---
  if (is.null(D_prior)) D_prior <- diag(p)
  if (is.null(Theta_init)) Theta <- matrix(0, K, K) else Theta <- as.matrix(Theta_init)
  # nu = -1.386 corresponds to q = logit^{-1}(-1.386) = 0.20, the Beta(1,4) prior mean
  if (is.null(nu_init)) nu <- matrix(-1, p, p) else nu <- as.matrix(nu_init)
  if (is.null(C_init)) C_arr <- array(diag(p), dim = c(p, p, K)) else C_arr <- C_init

  if (!all(dim(D_prior) == c(p, p))) stop("D_prior must be p x p.")
  if (!all(dim(Theta) == c(K, K))) stop("Theta_init must be K x K.")
  if (!all(dim(nu) == c(p, p))) stop("nu_init must be p x p.")
  if (!all(dim(C_arr) == c(p, p, K))) stop("C_init must be p x p x K.")

  # Ensure symmetry
  nu <- (nu + t(nu)) / 2
  diag(nu) <- 0
  Theta <- (Theta + t(Theta)) / 2
  diag(Theta) <- 0

  # --- Derived quantities ---
  b_prior_vec <- rep_len(b_prior, K)
  D_prior_arr <- array(D_prior, dim = c(p, p, K))
  S_arr <- array(NA_real_, dim = c(p, p, K))
  for (k in seq_len(K)) S_arr[, , k] <- S_list[[k]]

  n_vec <- as.numeric(n_vec)
  b_post <- b_prior_vec + n_vec
  D_post_arr <- array(NA_real_, dim = c(p, p, K))
  for (k in seq_len(K)) D_post_arr[, , k] <- D_prior_arr[, , k] + S_arr[, , k]

  # Initial adjacency from precision matrix
  adj_arr <- array(as.integer(abs(C_arr) > 1e-5), dim = c(p, p, K))
  for (k in seq_len(K)) diag(adj_arr[, , k]) <- 0L

  # --- Allocate storage ---
  nmc <- as.integer(nsave) * as.integer(thin)
  niter <- as.integer(burnin) + nmc

  C_save     <- array(0, dim = c(p, p, K, nsave))
  Sig_save   <- array(0, dim = c(p, p, K, nsave))
  adj_save   <- array(0L, dim = c(p, p, K, nsave))
  Theta_save <- array(0, dim = c(K, K, nsave))
  nu_save    <- array(0, dim = c(p, p, nsave))

  # Acceptance tracking
  ar_gamma       <- matrix(0, K, K)   # between-model (spike/slab toggle)
  ar_theta       <- matrix(0, K, K)   # within-model (slab move)
  n_within_model <- matrix(0, K, K)
  ar_nu          <- matrix(0, p, p)

  save_idx <- 0L

  # =========================================================================
  # MCMC main loop
  # =========================================================================
  for (iter in seq_len(niter)) {
    if (isTRUE(disp) && (iter %% 500L == 0L)) message("iter = ", iter)

    # -----------------------------------------------------------------------
    # Step 1: Update (G_k, Omega_k) for each group k  [Appendix A.1]
    # Uses Wang & Li (2012) exchange algorithm
    # -----------------------------------------------------------------------
    for (cur_graph in seq_len(K)) {
      Ck   <- C_arr[, , cur_graph]
      adjk <- adj_arr[, , cur_graph]

      # Sanitize adjacency
      adjk[is.na(adjk)] <- 0L
      storage.mode(adjk) <- "integer"
      diag(adjk) <- 0L
      adjk <- ((adjk + t(adjk)) > 0L) * 1L

      for (ii in 1:(p - 1)) {
        for (jj in (ii + 1):p) {

          # Step 1(a): Bernoulli proposal from MRF conditional (eq A.3)
          logH <- log_H_cpp(b_prior_vec[cur_graph],
                            D_prior_arr[, , cur_graph],
                            n_vec[cur_graph],
                            S_arr[, , cur_graph],
                            Ck, ii, jj)

          # MRF conditional: nu_ij + 2 * sum_{m != k} theta_{km} * g_{m,ij}
          # (Theta diagonal is 0, so self-term vanishes automatically)
          sim_term <- nu[ii, jj] + 2 * sum(Theta[cur_graph, ] * as.numeric(adj_arr[ii, jj, ]))

          w_logodds <- logH - sim_term

          if (is.na(w_logodds)) {
            next
          } else if (is.infinite(w_logodds)) {
            w <- if (w_logodds > 0) 0 else 1
          } else {
            # Numerically stable logistic
            if (w_logodds > 0) {
              ew <- exp(-w_logodds)
              w <- ew / (1 + ew)
            } else {
              ew <- exp(w_logodds)
              w <- 1 / (1 + ew)
            }
          }

          if (!is.finite(w) || is.na(w)) next

          current_ij  <- adjk[ii, jj]
          propose_ij  <- as.integer(stats::runif(1) < w)

          if (is.na(propose_ij) || is.na(current_ij)) next

          if (isTRUE(propose_ij != current_ij)) {
            # Step 1(b): Sample proposal C under prior via Wang-Li
            out_prior <- GWishart_NOij_Gibbs_cpp(
              b_prior_vec[cur_graph],
              D_prior_arr[, , cur_graph],
              adjk * 1L, Ck, ii, jj,
              as.integer(propose_ij),
              0L, 1L
            )
            if (isFALSE(out_prior$ok)) next
            C_prop <- out_prior$C

            # Step 2(b): Exchange algorithm acceptance ratio
            r2 <- log_GWishart_NOij_pdf_cpp(b_prior_vec[cur_graph],
                                             D_prior_arr[, , cur_graph],
                                             C_prop, ii, jj,
                                             as.integer(current_ij)) -
              log_GWishart_NOij_pdf_cpp(b_prior_vec[cur_graph],
                                        D_prior_arr[, , cur_graph],
                                        C_prop, ii, jj,
                                        as.integer(propose_ij))

            if (log(stats::runif(1)) < r2) {
              adjk[ii, jj] <- propose_ij
              adjk[jj, ii] <- propose_ij
              current_ij <- propose_ij
            }
          }

          # Step 2(c): Update C given current graph under posterior
          out_g <- GWishart_NOij_Gibbs_cpp(b_post[cur_graph],
                                            D_post_arr[, , cur_graph],
                                            adjk * 1L, Ck, ii, jj,
                                            as.integer(current_ij),
                                            0L, 0L)
          if (isFALSE(out_g$ok)) next
          Ck   <- out_g$C
          adjk <- out_g$adj
        }
      }

      # Step 3: BIPS clique update for full precision matrix
      out_bips <- GWishart_BIPS_maximumClique_cpp(b_post[cur_graph],
                                                    D_post_arr[, , cur_graph],
                                                    adjk * 1L, Ck, 0L, 1L)
      if (isFALSE(out_bips$ok)) next

      C_arr[, , cur_graph]   <- out_bips$C
      adj_arr[, , cur_graph] <- adjk
    }

    # -----------------------------------------------------------------------
    # Step 2: Update (theta_km, gamma_km)  [Appendix A.2]
    # Between-model (spike/slab toggle) + within-model (slab) moves
    # -----------------------------------------------------------------------
    for (k in 1:(K - 1)) {
      for (m in (k + 1):K) {

        # --- Between-model move ---
        theta_curr <- Theta[k, m]
        # If theta=0 (spike), propose from Gamma slab; if theta>0 (slab), propose 0
        theta_prop <- if (theta_curr == 0) {
          stats::rgamma(1, shape = alpha_prop, rate = beta_prop)
        } else {
          0
        }

        Theta_prop <- Theta
        Theta_prop[k, m] <- theta_prop
        Theta_prop[m, k] <- theta_prop

        # Sum over all edges: logC(Theta_curr, nu_ij) - logC(Theta_prop, nu_ij) +
        #   2*(theta_prop - theta_curr) * g_{k,ij} * g_{m,ij}
        sum_over_edges <- 0
        for (ii in 1:(p - 1)) {
          for (jj in (ii + 1):p) {
            sum_over_edges <- sum_over_edges +
              calc_mrf_logC(Theta,      nu[ii, jj]) +
              2 * (theta_prop - theta_curr) *
              as.numeric(adj_arr[ii, jj, k]) * as.numeric(adj_arr[ii, jj, m]) -
              calc_mrf_logC(Theta_prop, nu[ii, jj])
          }
        }

        # Log acceptance ratio (eqs A.5, A.6)
        # Prior: spike-and-slab mixture with Gamma(alpha, rate=beta) slab
        log_prior_prop <- if (theta_prop == 0) {
          log(1 - my_w)
        } else {
          log(my_w) + stats::dgamma(theta_prop, shape = alpha, rate = beta, log = TRUE)
        }
        log_prior_curr <- if (theta_curr == 0) {
          log(1 - my_w)
        } else {
          log(my_w) + stats::dgamma(theta_curr, shape = alpha, rate = beta, log = TRUE)
        }

        # Proposal kernel: Gamma(alpha_prop, rate=beta_prop)
        log_q_prop_given_curr <- if (theta_curr == 0) {
          stats::dgamma(theta_prop, shape = alpha_prop, rate = beta_prop, log = TRUE)
        } else {
          0
        }
        log_q_curr_given_prop <- if (theta_prop == 0) {
          stats::dgamma(theta_curr, shape = alpha_prop, rate = beta_prop, log = TRUE)
        } else {
          0
        }

        log_ar <- (log_prior_prop - log_prior_curr) + sum_over_edges +
                  (log_q_curr_given_prop - log_q_prop_given_curr)

        if (is.finite(log_ar) && log_ar > log(stats::runif(1))) {
          Theta[k, m] <- theta_prop
          Theta[m, k] <- theta_prop
          ar_gamma[k, m] <- ar_gamma[k, m] + 1 / niter
        }

        # --- Within-model (slab) move (eq A.7) ---
        if (Theta[k, m] != 0) {
          theta_curr2  <- Theta[k, m]
          theta_prop2  <- stats::rgamma(1, shape = alpha_prop, rate = beta_prop)
          n_within_model[k, m] <- n_within_model[k, m] + 1

          Theta_prop2 <- Theta
          Theta_prop2[k, m] <- theta_prop2
          Theta_prop2[m, k] <- theta_prop2

          sum_over_edges2 <- 0
          for (ii in 1:(p - 1)) {
            for (jj in (ii + 1):p) {
              sum_over_edges2 <- sum_over_edges2 +
                calc_mrf_logC(Theta,       nu[ii, jj]) +
                2 * (theta_prop2 - theta_curr2) *
                as.numeric(adj_arr[ii, jj, k]) * as.numeric(adj_arr[ii, jj, m]) -
                calc_mrf_logC(Theta_prop2, nu[ii, jj])
            }
          }

          # Target ratio: Gamma(alpha, rate=beta) prior
          log_target_ratio <- (stats::dgamma(theta_prop2, shape = alpha, rate = beta, log = TRUE) -
                                 stats::dgamma(theta_curr2, shape = alpha, rate = beta, log = TRUE)) +
            sum_over_edges2

          # Hastings correction for independence Gamma(alpha_prop, rate=beta_prop) proposal
          log_hastings <- stats::dgamma(theta_curr2, shape = alpha_prop, rate = beta_prop, log = TRUE) -
            stats::dgamma(theta_prop2, shape = alpha_prop, rate = beta_prop, log = TRUE)

          log_theta_ar <- log_target_ratio + log_hastings

          if (is.finite(log_theta_ar) && log_theta_ar > log(stats::runif(1))) {
            Theta[k, m] <- theta_prop2
            Theta[m, k] <- theta_prop2
            ar_theta[k, m] <- ar_theta[k, m] + 1
          }
        }
      }
    }

    # -----------------------------------------------------------------------
    # Step 3: Update nu_ij  [Appendix A.3]
    # Independent proposals q* ~ Beta(a_prop, b_prop), nu* = logit(q*)
    # -----------------------------------------------------------------------
    for (ii in 1:(p - 1)) {
      for (jj in (ii + 1):p) {
        q <- stats::rbeta(1, a_prop, b_prop)
        nu_prop <- log(q) - log(1 - q)

        # Log MH ratio (eq A.10)
        log_nu_ar <- (nu_prop - nu[ii, jj]) * (sum(adj_arr[ii, jj, ]) + a - a_prop) -
          (a + b - a_prop - b_prop) * log(1 + exp(nu_prop)) -
          calc_mrf_logC(Theta, nu_prop) +
          (a + b - a_prop - b_prop) * log(1 + exp(nu[ii, jj])) +
          calc_mrf_logC(Theta, nu[ii, jj])

        if (log_nu_ar > log(stats::runif(1))) {
          nu[ii, jj] <- nu_prop
          nu[jj, ii] <- nu_prop
          ar_nu[ii, jj] <- ar_nu[ii, jj] + 1 / niter
        }
      }
    }

    # -----------------------------------------------------------------------
    # Store posterior samples
    # -----------------------------------------------------------------------
    if (iter > burnin && ((iter - burnin) %% thin == 0L)) {
      save_idx <- save_idx + 1L
      if (save_idx <= nsave) {
        C_save[, , , save_idx]   <- C_arr
        for (k in seq_len(K)) {
          Sig_save[, , k, save_idx] <- solve(C_arr[, , k])
        }
        adj_save[, , , save_idx] <- adj_arr * 1L
        Theta_save[, , save_idx] <- Theta
        nu_save[, , save_idx]    <- nu
      }
    }
  }

  # --- Post-process acceptance rates ---
  for (k in 1:(K - 1)) {
    for (m in (k + 1):K) {
      if (n_within_model[k, m] > 0) {
        ar_theta[k, m] <- ar_theta[k, m] / n_within_model[k, m]
      }
    }
  }

  structure(list(
    K = K,
    p = p,
    C_save     = C_save,
    Sig_save   = Sig_save,
    adj_save   = adj_save,
    Theta_save = Theta_save,
    nu_save    = nu_save,
    ar_gamma   = ar_gamma,
    ar_theta   = ar_theta,
    ar_nu      = ar_nu,
    hyper = list(alpha = alpha, beta = beta, a = a, b = b, w = my_w,
                 alpha_prop = alpha_prop, beta_prop = beta_prop,
                 a_prop = a_prop, b_prop = b_prop),
    call = match.call()
  ), class = "multiggm_fit")
}
