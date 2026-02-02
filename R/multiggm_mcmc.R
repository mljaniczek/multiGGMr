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

  #number of sample groups
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

#' Internal single-chain engine
#'
#' Implements the MATLAB routine \code{MCMC_multiple_graphs.m} using the Wang–Li
#' G-Wishart kernels (ported to RcppArmadillo).
#'
#' @note internal
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
        w = 0.1,
        alpha_prop = 1, beta_prop = 1
      )
    }

    # Then refer to hyper$... throughout:
    b_prior <- hyper$b_prior
    D_prior <- hyper$D_prior
    a <- hyper$a; b <- hyper$b
    alpha <- hyper$alpha; beta <- hyper$beta
    my_w <- hyper$w
    alpha_prop <- hyper$alpha_prop
    beta_prop <- hyper$beta_prop

  stopifnot(is.list(S_list), length(S_list) >= 1L)
  K <- length(S_list)
  p <- nrow(S_list[[1L]])
  if (is.null(D_prior)) D_prior <- diag(p)
  if (is.null(Theta_init)) Theta <- matrix(0, K, K) else Theta <- as.matrix(Theta_init)
  # why is nu starting out as all -1?
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
  # TODO should these be a and be or alpha and beta? or static within?
  a_prop <- 2
  b_prop <- 4

  nmc <- as.integer(nsave) * as.integer(thin)
  niter <- as.integer(burnin) + nmc

  # set up matrices for return values
  C_save <- array(0, dim = c(p, p, K, nsave))
  Sig_save <- array(0, dim = c(p, p, K, nsave))
  adj_save <- array(0L, dim = c(p, p, K, nsave))
  Theta_save <- array(0, dim = c(K, K, nsave))
  nu_save <- array(0, dim = c(p, p, nsave))

  # acceptance rate for gamma based on number of between model moves
  ar_gamma <- matrix(0, K, K)

  # acceptance rate for theta based on number of within model moves
  ar_theta <- matrix(0, K, K)
  n_within_model <- matrix(0, K, K)
  ar_nu <- matrix(0, p, p)

  # --- bookkeeping for spike/slab (between-model) toggles ---
  toggle_attempt <- matrix(0L, K, K)
  toggle_accept  <- matrix(0L, K, K)

  # bookkeeping for within-slab Theta proposals (optional, but useful)
  theta_attempt  <- matrix(0L, K, K)
  theta_accept   <- matrix(0L, K, K)
  # DEBUG: theta toggle diagnostics (do not leave on forever)
  dbg_theta_prop_pos <- 0L      # number of times theta_prop > 0 is proposed
  dbg_theta_logar_finite <- 0L  # number of times log_ar is finite
  dbg_theta_logar_gt0 <- 0L     # number of times log_ar > 0 (would accept with prob > 0.5)
  dbg_theta_accept <- 0L
  dbg_toggle_up_attempt <- 0L   # theta_curr==0, propose >0
  dbg_toggle_up_accept  <- 0L
  dbg_toggle_dn_attempt <- 0L   # theta_curr>0, propose 0
  dbg_toggle_dn_accept  <- 0L
  dbg_logar_up_gt0 <- 0L
  dbg_logar_dn_gt0 <- 0L
  dbg_w <- 0L
  dbg_r2_infinite <- 0L


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

      # sample off-diagonal elements for each graph
      for (ii in 1:(p - 1)) {
        for (jj in (ii + 1):p) {
      #Step 1(a) Bernoulli proposal
          # log odds (MATLAB uses log_H - (nu + 2 Theta*adj))
          logH <- log_H_cpp(b_prior_vec[cur_graph],
                            D_prior_arr[, , cur_graph],
                            n_vec[cur_graph],
                            S_arr[, , cur_graph],
                            Ck, ii, jj)

          # similarity penalty term
          #?? where dis? is this the bernoulli prob?
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

          # indicates whether edge [[i,j]] is in G
          current_ij <- adjk[ii, jj]

          # indicates whether edge [[i,j]] is in G'
          # proposal will be 1 if runif(1)<w (i.e. will be an edge) with probability w
          propose_ij <- as.integer(stats::runif(1) < w)  # force 0/1 integer

          ## TODO add check here to collect the ws to see if this is a problem?
          dbg_w <- dbg_w + as.integer(stats::runif(1) < w)

          if (is.na(propose_ij) || is.na(current_ij)) next
          if (isTRUE(propose_ij != current_ij)) {
            # Step 1(b): sample proposal C under prior
            ## function adapted from from Wang/Li paper and is in wangli.cpp
            out_prior <- GWishart_NOij_Gibbs_cpp(
              b_prior_vec[cur_graph],
              D_prior_arr[, , cur_graph],
              adjk * 1L, Ck, ii, jj,
              as.integer(propose_ij),
              0L, 1L
            )
            # soft move on if it didn't work in one iter
            if (isFALSE(out_prior$ok)) next
            # otherwise collect it and move on
            C_prop <- out_prior$C

            # Step 2(b) from expression at top of p.189 of Wang/Li paper (function adapted from them, contained in wangli.cpp)
            ### TODO this is infinite!! problem??
            # if it's alwasy infinite then  acceptance will always happen
            r2 <- log_GWishart_NOij_pdf_cpp(b_prior_vec[cur_graph],
                                           D_prior_arr[, , cur_graph],
                                           C_prop, ii, jj,
                                           as.integer(current_ij)) -
              log_GWishart_NOij_pdf_cpp(b_prior_vec[cur_graph],
                                        D_prior_arr[, , cur_graph],
                                        C_prop, ii, jj,
                                        as.integer(propose_ij))

            # acceptance rate alpha = min(1, e^r2)
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
            # reject/skip update if it doesn't work in one iter rather than burn everything down
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
        # Numerical failure inside BIPS update; reject/skip and keep current state instead of burning everything downnnn
        next
      }
      C_arr[, , cur_graph] <- out_bips$C
      Sigk <- out_bips$Sig
      adj_arr[, , cur_graph] <- adjk

      if (iter > burnin && ((iter - burnin) %% thin == 0L)) {
        # Sig stored with thinning
        # will be copied below at save time
        # MJ TODO:: Double check this happens?? why is this blank
      }
    }

    # --- update Theta (network relatedness): MATLAB-faithful Gamma proposals ---
    for (k in 1:(K - 1)) {
      for (m in (k + 1):K) {
        # -----------------------------
        # Between-model (spike/slab) move
        # MATLAB: if theta==0 propose Gamma; else propose 0
        # -----------------------------
        theta_curr <- Theta[k, m]
        # DEBUG
        if (theta_curr == 0) {
          dbg_toggle_up_attempt <- dbg_toggle_up_attempt + 1L
        } else {
          dbg_toggle_dn_attempt <- dbg_toggle_dn_attempt + 1L
        }
        # end debug
        theta_prop <- if (theta_curr == 0) stats::rgamma(1, shape = alpha_prop, scale = beta_prop) else 0

        # bookkeeping: attempted toggle
        toggle_attempt[k, m] <- toggle_attempt[k, m] + 1L
        if (theta_prop > 0) dbg_theta_prop_pos <- dbg_theta_prop_pos + 1L

        Theta_prop <- Theta
        Theta_prop[k, m] <- theta_prop
        Theta_prop[m, k] <- theta_prop

        # get terms that are a sum over all edges on log scale
        sum_over_edges <- 0
        for (ii in 1:(p - 1)) {
          for (jj in (ii + 1):p) {
            sum_over_edges <- sum_over_edges +
              ### is this where there's a discrepancy??
              # TODO double check calc_mrf_logC function
              calc_mrf_logC(Theta,      nu[ii, jj]) +
              2 * (theta_prop - theta_curr) *
              as.numeric(adj_arr[ii, jj, k]) * as.numeric(adj_arr[ii, jj, m]) -
              calc_mrf_logC(Theta_prop, nu[ii, jj])
          }
        } # in my testing one iteration this val was 183.01... is that too big?

        # MATLAB-style log acceptance ratio



        ## MJ modifying original to more closely follow matlab code

        # if (theta_prop ==0) {
        #   log_ar = alpha_prop*log(beta_prop) - log(gamma(alpha_prop)) +
        #     log(gamma(alpha)) - alpha*log(beta) - (alpha - alpha_prop)*log(theta_curr) +
        #     sum_over_edges + log(1-my_w) - log(my_w)
        # } else {
        #   log_ar = alpha*log(beta) - log(gamma(alpha)) + log(gamma(alpha_prop)) -
        #     alpha_prop*log(beta_prop) - (alpha-alpha_prop)*log(theta_prop)-
        #     (beta-beta_prop)*theta_prop + sum_over_edges + log(my_w) - log(1-my_w)
        # } ## then the log_ar was 178.6 mainly bc of the sum of edges

        #welp the above didn't work so going back to other
        # Mixture prior on theta:
        log_prior_prop <- if (theta_prop == 0) log(1 - my_w) else log(my_w) + dgamma(theta_prop, shape = alpha, scale = beta, log = TRUE)
        log_prior_curr <- if (theta_curr == 0) log(1 - my_w) else log(my_w) + dgamma(theta_curr, shape = alpha, scale = beta, log = TRUE)

        # Toggle proposal kernel:
        log_q_prop_given_curr <- if (theta_curr == 0) dgamma(theta_prop, shape = alpha_prop, scale = beta_prop, log = TRUE) else 0
        log_q_curr_given_prop <- if (theta_prop == 0) dgamma(theta_curr, shape = alpha_prop, scale = beta_prop, log = TRUE) else 0

        log_ar <- (log_prior_prop - log_prior_curr) + sum_over_edges + (log_q_curr_given_prop - log_q_prop_given_curr)

       #  # Prior mixture on theta (spike at 0 with prob 1-w; slab Gamma(alpha,beta))
       #  log_prior_prop <- if (theta_prop == 0) log(1 - my_w) else (log(my_w) + stats::dgamma(theta_prop, shape = alpha, scale = beta, log = TRUE))
       #  ## another fix possibly - this previously was both log(1-my_w)
       #  log_prior_curr <- if (theta_curr == 0) log(1- my_w) else (log(my_w) + stats::dgamma(theta_curr, shape = alpha, scale = beta, log = TRUE))
       #
       #  # Proposal densities: q(theta_prop | theta_curr)
       #  # if leaving spike: propose from Gamma(alpha_prop,beta_prop)
       #  # if leaving slab: propose 0 deterministically
       #  log_q_prop_given_curr <- if (theta_curr == 0) stats::dgamma(theta_prop, shape = alpha_prop, scale = beta_prop, log = TRUE) else 0
       #  log_q_curr_given_prop <- if (theta_prop == 0) stats::dgamma(theta_curr, shape = alpha_prop, scale = beta_prop, log = TRUE) else 0
       #
       #  log_ar <- (log_prior_prop - log_prior_curr) + sum_over_edges + (log_q_curr_given_prop - log_q_prop_given_curr)
       # #debug below
         if (is.finite(log_ar)) dbg_theta_logar_finite <- dbg_theta_logar_finite + 1L
        if (is.finite(log_ar) && log_ar > 0) dbg_theta_logar_gt0 <- dbg_theta_logar_gt0 + 1L
        if (theta_curr == 0 && is.finite(log_ar) && log_ar > 0) dbg_logar_up_gt0 <- dbg_logar_up_gt0 + 1L
        if (theta_curr != 0 && is.finite(log_ar) && log_ar > 0) dbg_logar_dn_gt0 <- dbg_logar_dn_gt0 + 1L
#end debug

        # accept proposal with given prob
        # will this always be accepted bc of log_ar bweing so big?
        if (is.finite(log_ar) && log_ar > log(stats::runif(1))) {
          Theta[k, m] <- theta_prop
          Theta[m, k] <- theta_prop

          # bookkeeping: toggle accepted
          toggle_accept[k, m] <- toggle_accept[k, m] + 1L
          # increment acceptance scale
          ar_gamma[k, m] <- ar_gamma[k, m] + 1 / niter
          dbg_theta_accept <- dbg_theta_accept + 1L
          if (theta_curr == 0) dbg_toggle_up_accept <- dbg_toggle_up_accept + 1L
          if (theta_curr != 0) dbg_toggle_dn_accept <- dbg_toggle_dn_accept + 1L
          #end bookkeeping
        }

        # -----------------------------
        # Within-model (slab) move
        # MATLAB: if theta!=0 propose Gamma(alpha_prop,beta_prop) (independence proposal)
        # -----------------------------
        if (Theta[k, m] != 0) {
          # TODO find the n_within_model(k,m)

          theta_curr2 <- Theta[k, m]
          # so this is basically saying, if Theta is not zero, sample from gamma distribution.
          # average rgamma() when alpha_prop = beta_prop = 1 should be 1
          theta_prop2 <- stats::rgamma(1, shape = alpha_prop, scale = beta_prop)

          # bookkeeping: within-slab attempt
          theta_attempt[k, m] <- theta_attempt[k, m] + 1L
          n_within_model[k, m] <- n_within_model[k, m] + 1

          Theta_prop2 <- Theta
          Theta_prop2[k, m] <- theta_prop2
          Theta_prop2[m, k] <- theta_prop2

          # get terms that are a sum over all edges on log scale
          sum_over_edges2 <- 0
          for (ii in 1:(p - 1)) {
            for (jj in (ii + 1):p) {
              sum_over_edges2 <- sum_over_edges2 +
                # MJ TODO double check this function - I think this is the root of the problem?
                calc_mrf_logC(Theta,       nu[ii, jj]) +
                2 * (theta_prop2 - theta_curr2) *
                as.numeric(adj_arr[ii, jj, k]) * as.numeric(adj_arr[ii, jj, m]) -
                calc_mrf_logC(Theta_prop2, nu[ii, jj])
            }
          }

          # Target ratio in slab: Gamma(alpha,beta) prior for theta * exp(sum_over_edges2)
          # log_target_ratio <- (stats::dgamma(theta_prop2, shape = alpha, scale = beta, log = TRUE) -
          #                        stats::dgamma(theta_curr2, shape = alpha, scale = beta, log = TRUE)) +
          #   sum_over_edges2
          #
          # # Hastings for independence proposal q(.) = Gamma(alpha_prop,beta_prop)
          # log_hastings <- stats::dgamma(theta_curr2, shape = alpha_prop, scale = beta_prop, log = TRUE) -
          #   stats::dgamma(theta_prop2, shape = alpha_prop, scale = beta_prop, log = TRUE)
          #
          # log_theta_ar <- log_target_ratio + log_hastings

          #MJ redoing the above to match the matlab code better
          # log_theta_ar = (alpha - alpha_prop)*(log(theta_prop2) - log(theta_curr2)) +
          #   (beta - beta_prop)*(theta_curr2 - theta_prop2) + sum_over_edges2

          # once again my try didnt work so going back to try below instead:
          log_target_ratio <- (dgamma(theta_prop2, shape=alpha, scale=beta, log=TRUE) -
                                 dgamma(theta_curr2, shape=alpha, scale=beta, log=TRUE)) +
            sum_over_edges2

          log_hastings <- dgamma(theta_curr2, shape=alpha_prop, scale=beta_prop, log=TRUE) -
            dgamma(theta_prop2, shape=alpha_prop, scale=beta_prop, log=TRUE)

          log_theta_ar <- log_target_ratio + log_hastings


          # accept proposal with given probability
          if (is.finite(log_theta_ar) && log_theta_ar > log(stats::runif(1))) {
            Theta[k, m] <- theta_prop2
            Theta[m, k] <- theta_prop2
            # track number of proposals accepted
            ar_theta[k, m] <- ar_theta[k, m] + 1
  #TODO go through and remove redundant bookkeeping
            # bookkeeping: within-slab accept
            theta_accept[k, m] <- theta_accept[k, m] + 1L
          }
        }

      }
    }
    pos <- 0; tot <- 0 # TODO label these

    # generate independent proposals for q from beta(aprop, bprop) density
    # --- update nu ---
    for (ii in 1:(p - 1)) {
      for (jj in (ii + 1):p) {
        q <- stats::rbeta(1, a_prop, b_prop)
        nu_prop <- log(q) - log(1 - q)

        # calculate MH ratio on log scale

        log_nu_ar <- (nu_prop - nu[ii, jj]) * (sum(adj_arr[ii, jj, ]) + a - a_prop) -
          (a + b - a_prop - b_prop) * log(1 + exp(nu_prop)) -
          # TODO check this function
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
      # number of accepted toggles (should match toggle_accept sum)
    dbg_theta = list(
      theta_prop_pos = dbg_theta_prop_pos,
      logar_finite   = dbg_theta_logar_finite,
      logar_gt0      = dbg_theta_logar_gt0,
      accept         = dbg_theta_accept,
      dbg_toggle_up_attempt = dbg_toggle_up_attempt,   # theta_curr==0, propose >0
      dbg_toggle_up_accept  =dbg_toggle_up_accept,
      dbg_toggle_dn_attempt =dbg_toggle_dn_attempt,   # theta_curr>0, propose 0
      dbg_toggle_dn_accept  =dbg_toggle_dn_accept,
      dbg_logar_up_gt0 =dbg_logar_up_gt0,
      dbg_logar_dn_gt0 =dbg_logar_dn_gt0,
      dbg_w = dbg_w
    ),
    hyper = list(alpha = alpha, beta = beta, a = a, b = b, my_w = my_w,
                 alpha_prop = alpha_prop, beta_prop = beta_prop, a_prop = a_prop, b_prop = b_prop),
    call = match.call()
  ), class = "multiggm_fit")
}

