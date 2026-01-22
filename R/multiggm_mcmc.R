#' Fit multi-graph Gaussian graphical model via MCMC
#'
#' This is the package-level entry point that mirrors the MATLAB workflow
#' (MCMC_multiple_graphs.m) while providing an R-native interface and
#' optional parallel execution for multiple chains.
#'
#' @param S_list List of K sample covariance matrices (each p x p).
#' @param n_vec Integer vector of length K with sample sizes.
#' @param burnin Integer burn-in iterations.
#' @param nsave Integer number of saved posterior draws per chain (after thinning).
#' @param thin Integer thinning interval.
#' @param nchains Number of independent MCMC chains.
#' @param parallel Logical; if TRUE, run chains in parallel using \pkg{parallel}.
#' @param ncores Number of cores/workers to use when \code{parallel=TRUE}.
#' @param seed Optional integer seed. If provided, chains use deterministic offsets.
#' @param ... Additional arguments forwarded to the single-chain engine.
#'
#' @return An object of class \code{"multiggm_fit"} (single chain) or
#'   \code{"multiggm_fit_list"} (multiple chains).
#' @export
multiggm_mcmc <- function(S_list, n_vec,
                          burnin = 5e3, nsave = 1e3, thin = 1,
                          nchains = 1,
                          parallel = FALSE, ncores = max(1L, parallel::detectCores() - 1L),
                          seed = NULL,
                          ...) {

  if (!is.list(S_list) || length(S_list) < 1L) stop("S_list must be a non-empty list of covariance matrices.")
  K <- length(S_list)
  n_vec <- as.integer(n_vec)
  if (length(n_vec) != K) stop("n_vec must have length equal to length(S_list).")
  if (!is.null(seed)) seed <- as.integer(seed)

  args <- list(...)
  chain_fun <- function(chain_id) {
    if (!is.null(seed)) set.seed(seed + 1000L * (chain_id - 1L))
    .multiggm_mcmc_single(S_list = S_list, n_vec = n_vec,
                         burnin = burnin, nsave = nsave, thin = thin,
                         chain_id = chain_id, ...)
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
      # Export needed objects/functions to workers
      parallel::clusterExport(cl, varlist = c("S_list", "n_vec", "burnin", "nsave", "thin", "seed"),
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

#' Internal single-chain engine (to be implemented)
#'
#' This placeholder exists to lock the API and the parallel orchestration.
#' The full Wang-Li G-Wishart kernels will be wired here.
#'
#' @keywords internal

#' Internal single-chain engine
#'
#' Implements the MATLAB routine \code{MCMC_multiple_graphs.m} using the Wang–Li
#' G-Wishart kernels (ported to RcppArmadillo).
#'
#' @keywords internal
.multiggm_mcmc_single <- function(S_list, n_vec, burnin, nsave, thin,
                                 chain_id = 1L,
                                 Theta_init = NULL,
                                 nu_init = NULL,
                                 C_init = NULL,
                                 disp = FALSE,
                                 hyper = NULL,
                                seed = NULL, ...) {

    # Default hyperparameters if not provided
    if (is.null(hyper)) {
      p <- ncol(S_list[[1]])
      hyper <- list(
        b_prior = 3,
        D_prior = diag(p),
        a = 1, b = 1,
        alpha = 2, beta = 2,
        w = 0.1
      )
    }

    # Then refer to hyper$... throughout:
    b_prior <- hyper$b_prior
    D_prior <- hyper$D_prior
    a <- hyper$a; b <- hyper$b
    alpha <- hyper$alpha; beta <- hyper$beta
    my_w <- hyper$w

  stopifnot(is.list(S_list), length(S_list) >= 1L)
  K <- length(S_list)
  p <- nrow(S_list[[1L]])
  if (is.null(D_prior)) D_prior <- diag(p)
  if (is.null(Theta_init)) Theta <- matrix(0, K, K) else Theta <- as.matrix(Theta_init)
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

  # Priors per group
  b_prior_vec <- rep_len(b_prior, K)
  b_prior_mat <- matrix(b_prior_vec, nrow = 1L) # for MATLAB-compat shapes
  D_prior_arr <- array(D_prior, dim = c(p, p, K))
  S_arr <- array(NA_real_, dim = c(p, p, K))
  for (k in seq_len(K)) S_arr[, , k] <- S_list[[k]]

  n_vec <- as.numeric(n_vec)
  b_post <- b_prior_vec + n_vec
  D_post_arr <- array(NA_real_, dim = c(p, p, K))
  for (k in seq_len(K)) D_post_arr[, , k] <- D_prior_arr[, , k] + S_arr[, , k]

  # Initial adjacency inferred from C
  adj_arr <- array(as.integer(abs(C_arr) > 1e-5), dim = c(p, p, K))
  for (k in seq_len(K)) diag(adj_arr[, , k]) <- 0L


  # Proposal hyperparameters (as in MATLAB)
  alpha_prop <- 1
  beta_prop <- 1
  a_prop <- 2
  b_prop <- 4

  nmc <- as.integer(nsave) * as.integer(thin)
  niter <- as.integer(burnin) + nmc

  C_save <- array(0, dim = c(p, p, K, nsave))
  Sig_save <- array(0, dim = c(p, p, K, nsave))
  adj_save <- array(0L, dim = c(p, p, K, nsave))
  Theta_save <- array(0, dim = c(K, K, nsave))
  nu_save <- array(0, dim = c(p, p, nsave))

  ar_gamma <- matrix(0, K, K)
  ar_theta <- matrix(0, K, K)
  n_within_model <- matrix(0, K, K)
  ar_nu <- matrix(0, p, p)

  # --- bookkeeping for spike/slab (between-model) toggles ---
  toggle_attempt <- matrix(0L, K, K)
  toggle_accept  <- matrix(0L, K, K)

  # bookkeeping for within-slab Theta proposals (optional, but useful)
  theta_attempt  <- matrix(0L, K, K)
  theta_accept   <- matrix(0L, K, K)

  save_idx <- 0L

  for (iter in seq_len(niter)) {
    if (isTRUE(disp) && (iter %% 500L == 0L)) message("iter = ", iter)

    # --- update each group graph + precision ---
    for (cur_graph in seq_len(K)) {
      Ck <- C_arr[, , cur_graph]
      adjk <- adj_arr[, , cur_graph]

      # sanitize adjacency (0/1, symmetric, zero diagonal, no NA)
      adjk[is.na(adjk)] <- 0L
      storage.mode(adjk) <- "integer"
      diag(adjk) <- 0L
      adjk <- ((adjk + t(adjk)) > 0L) * 1L

      for (ii in 1:(p - 1)) {
        for (jj in (ii + 1):p) {

          # log odds (MATLAB uses log_H - (nu + 2 Theta*adj))
          logH <- log_H_cpp(b_prior_vec[cur_graph],
                            D_prior_arr[, , cur_graph],
                            n_vec[cur_graph],
                            S_arr[, , cur_graph],
                            Ck, ii, jj)

          # similarity penalty term
          sim_term <- nu[ii, jj] + 2 * sum(Theta[cur_graph, ] * as.numeric(adj_arr[ii, jj, ]))
          # compute w safely; if undefined, skip this edge update
          w_logodds <- logH - sim_term

          # Robust handling:
          # - NA -> skip this (ii,jj) update
          # - +Inf -> w = 0
          # - -Inf -> w = 1
          # - finite -> logistic transform
          # if (is.na(w_logodds)) {
          #   message("NA w_logodds at k=", cur_graph, " ii=", ii, " jj=", jj,
          #           " logH=", logH, " sim_term=", sim_term,
          #           " nu=", nu[ii,jj])
          #   next
          # }
          if (is.na(w_logodds)) {
            next
          } else if (is.infinite(w_logodds)) {
            w <- if (w_logodds > 0) 0 else 1
          } else {
            # numerically stable logistic
            if (w_logodds > 0) {
              ew <- exp(-w_logodds)
              w <- ew / (1 + ew)
            } else {
              ew <- exp(w_logodds)
              w <- 1 / (1 + ew)
            }
          }

          if (!is.finite(w) || is.na(w)) next

          current_ij <- adjk[ii, jj]
          propose_ij <- as.integer(stats::runif(1) < w)  # force 0/1 integer

          if (is.na(propose_ij) || is.na(current_ij)) next
          if (isTRUE(propose_ij != current_ij)) {
            # Step 1(b): sample proposal C under prior
            out_prior <- GWishart_NOij_Gibbs_cpp(
              b_prior_vec[cur_graph],
              D_prior_arr[, , cur_graph],
              adjk * 1L, Ck, ii, jj,
              as.integer(propose_ij),
              0L, 1L
            )
            if (isFALSE(out_prior$ok)) next
            C_prop <- out_prior$C

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

          # Step 2(c): update C given current graph under posterior
          out_g <- GWishart_NOij_Gibbs_cpp(b_post[cur_graph],
                                          D_post_arr[, , cur_graph],
                                          adjk*1L, Ck, ii, jj,
                                          as.integer(current_ij),
                                          0L, 0L)
          if (isFALSE(out_g$ok)) {
            # reject/skip update
            next
          }
          Ck <- out_g$C
          adjk <- out_g$adj
        }
      }

      # Step 3: update C and Sigma given graph
      out_bips <- GWishart_BIPS_maximumClique_cpp(b_post[cur_graph],
                                                  D_post_arr[, , cur_graph],
                                                  adjk*1L, Ck, 0L, 1L)
      if (isFALSE(out_bips$ok)) {
        # Numerical failure inside BIPS update; reject/skip and keep current state
        next
      }
      C_arr[, , cur_graph] <- out_bips$C
      Sigk <- out_bips$Sig
      adj_arr[, , cur_graph] <- adjk

      if (iter > burnin && ((iter - burnin) %% thin == 0L)) {
        # Sig stored with thinning
        # will be copied below at save time
      }
    }

    # --- update Theta (network relatedness) ---
    for (k in 1:(K - 1)) {
      for (m in (k + 1):K) {

        # Between-model (spike/slab) move
        theta_prop <- if (Theta[k, m] == 0) stats::rgamma(1, shape = alpha_prop, scale = beta_prop) else 0
        # bookkeeping
        toggle_attempt[k, m] <- toggle_attempt[k, m] + 1L

        Theta_prop <- Theta
        Theta_prop[k, m] <- theta_prop
        Theta_prop[m, k] <- theta_prop

        sum_over_edges <- 0
        for (ii in 1:(p - 1)) {
          for (jj in (ii + 1):p) {
            # log C(Theta, nu) terms
            sum_over_edges <- sum_over_edges +
              calc_mrf_logC(Theta, nu[ii, jj]) +
              2 * (theta_prop - Theta[m, k]) * as.numeric(adj_arr[ii, jj, k]) * as.numeric(adj_arr[ii, jj, m]) -
              calc_mrf_logC(Theta_prop, nu[ii, jj])
          }
        }

        if (theta_prop == 0) {
          log_ar <- alpha_prop * log(beta_prop) - lgamma(alpha_prop) +
            lgamma(alpha) - alpha * log(beta) -
            (alpha - alpha_prop) * log(Theta[m, k]) +
            (beta - beta_prop) * (Theta[m, k]) + sum_over_edges +
            log(1 - my_w) - log(my_w)
        } else {
          log_ar <- alpha * log(beta) - lgamma(alpha) +
            lgamma(alpha_prop) - alpha_prop * log(beta_prop) -
            (alpha - alpha_prop) * log(theta_prop) -
            (beta - beta_prop) * theta_prop + sum_over_edges +
            log(my_w) - log(1 - my_w)
        }

        if (log_ar > log(stats::runif(1))) {
          Theta[k, m] <- theta_prop
          Theta[m, k] <- theta_prop
          toggle_accept[k, m] <- toggle_accept[k, m] + 1L #bookkeeping
          ar_gamma[k, m] <- ar_gamma[k, m] + 1 / niter
        }

        # Within-model move (only if slab)
        if (Theta[k, m] != 0) {
          theta_attempt[k, m] <- theta_attempt[k, m] + 1L #bookeeping
          n_within_model[k, m] <- n_within_model[k, m] + 1
          # log-random-walk proposal (local move)
          sd_rw <- .6
          theta_curr <- Theta[k, m]
          theta_prop2 <- theta_curr * exp(stats::rnorm(1, mean = 0, sd = sd_rw))

          Theta_prop2 <- Theta
          Theta_prop2[k, m] <- theta_prop2
          Theta_prop2[m, k] <- theta_prop2

          sum_over_edges2 <- 0
          for (ii in 1:(p - 1)) {
            for (jj in (ii + 1):p) {
              sum_over_edges2 <- sum_over_edges2 +
                calc_mrf_logC(Theta, nu[ii, jj]) +
                2 * (theta_prop2 - theta_curr) *
                as.numeric(adj_arr[ii, jj, k]) * as.numeric(adj_arr[ii, jj, m]) -
                calc_mrf_logC(Theta_prop2, nu[ii, jj])
            }
          }

          # Prior ratio under Gamma(alpha, rate=beta):  (alpha-1)log(theta) - beta*theta
          log_prior_ratio <- (alpha - 1) * (log(theta_prop2) - log(theta_curr)) -
            beta * (theta_prop2 - theta_curr)

          # Hastings correction for log-normal RW: q(curr|prop)/q(prop|curr) = prop/curr
          log_hastings <- log(theta_prop2) - log(theta_curr)

          log_theta_ar <- log_prior_ratio + sum_over_edges2 + log_hastings

          if (log_theta_ar > log(stats::runif(1))) {
            Theta[k, m] <- theta_prop2
            Theta[m, k] <- theta_prop2
            ar_theta[k, m] <- ar_theta[k, m] + 1
            # bookeeping
            theta_accept[k, m] <- theta_accept[k, m] + 1L #bookeeping
          }
        }
      }
    }
    pos <- 0; tot <- 0
    # --- update nu ---
    for (ii in 1:(p - 1)) {
      for (jj in (ii + 1):p) {
        q <- stats::rbeta(1, a_prop, b_prop)
        nu_prop <- log(q) - log(1 - q)

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
        tot <- tot + 1
        if (log_nu_ar > 0) pos <- pos + 1
      }
    }
    message("nu: fraction log_nu_ar > 0 = ", pos/tot)
    # --- retain sample ---
    if (iter > burnin && ((iter - burnin) %% thin == 0L)) {
      save_idx <- save_idx + 1L
      if (save_idx <= nsave) {
        C_save[, , , save_idx] <- C_arr
        for (k in seq_len(K)) {
          Sig_save[, , k, save_idx] <- solve(C_arr[, , k])
        }
        adj_save[, , , save_idx] <- adj_arr * 1L
        Theta_save[, , save_idx] <- Theta
        nu_save[, , save_idx] <- nu
      }
    }
  }

  # Post-process acceptance rates
  for (k in 1:(K - 1)) {
    for (m in (k + 1):K) {
      if (n_within_model[k, m] > 0) ar_theta[k, m] <- ar_theta[k, m] / n_within_model[k, m]
    }
  }

  structure(list(
    K = K,
    p = p,
    C_save = C_save,
    Sig_save = Sig_save,
    adj_save = adj_save,
    Theta_save = Theta_save,
    nu_save = nu_save,
    ar_gamma = ar_gamma,
    ar_theta = ar_theta,
    ar_nu = ar_nu,
    toggle_attempt = toggle_attempt,
    toggle_accept  = toggle_accept,
    theta_attempt  = theta_attempt,
    theta_accept   = theta_accept,
    hyper = list(alpha = alpha, beta = beta, a = a, b = b, my_w = my_w,
                 alpha_prop = alpha_prop, beta_prop = beta_prop, a_prop = a_prop, b_prop = b_prop),
    call = match.call()
  ), class = "multiggm_fit")
}

