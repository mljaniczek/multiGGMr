
#' Compute posterior partial correlations from a multiggm_fit object
#'
#' Extracts the full array of posterior partial correlation draws from a
#' fitted multi-GGM model. For each saved MCMC iteration and group, the
#' precision matrix is converted to a partial correlation matrix.
#'
#' @param fit A \code{multiggm_fit} object (single chain) or
#'   \code{multiggm_fit_list} object (multiple chains).
#' @param use Character; which posterior array to compute partial correlations
#'   from:
#'   \describe{
#'     \item{\code{"C_save"}}{(default) Use the precision matrix draws
#'       directly. Faster and recommended.}
#'     \item{\code{"Sig_save"}}{Use the covariance matrix draws, inverting
#'       each to get precision first. Slower; mainly for cross-checking.}
#'   }
#' @param chain \code{NULL} (default) pools the draws of all chains of a
#'   \code{multiggm_fit_list}; an integer vector selects chains.
#'
#' @return A numeric array of dimension \code{[p, p, K, nsave]}. Entry
#'   \code{pcor[i, j, k, s]} is the partial correlation between variables
#'   \code{i} and \code{j} in group \code{k} at saved iteration \code{s}.
#'   Diagonal entries are 1. Off-diagonal entries are in \code{[-1, 1]}.
#'
#' @details
#' The partial correlation between variables \code{i} and \code{j} given
#' all other variables is defined as:
#' \deqn{\rho_{ij|rest} = -\Omega_{ij} / \sqrt{\Omega_{ii} \Omega_{jj}}}
#' where \eqn{\Omega} is the precision matrix.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' pcor_draws <- posterior_pcor(fit)
#' # Posterior mean partial correlation for group 1
#' pcor_mean <- apply(pcor_draws[,,1,], c(1,2), mean)
#'
#' @seealso [precision_to_pcor()], [posterior_ci()], [fitted.multiggm_fit()]
#' @export
posterior_pcor <- function(fit, use = c("C_save", "Sig_save"), chain = NULL) {
  use <- match.arg(use)
  arrs <- lapply(.fit_chains(fit, chain), .posterior_pcor_one, use = use)
  if (length(arrs) == 1L) return(arrs[[1]])
  d <- dim(arrs[[1]])
  array(unlist(arrs, use.names = FALSE),
        dim = c(d[1:3], sum(vapply(arrs, function(x) dim(x)[4], numeric(1)))))
}

.posterior_pcor_one <- function(fit, use) {
  if (is.null(fit[[use]]))
    stop("This fit has no ", use, " draws. mGHS fits store precision draws ",
         "only (C_save), and none when save_draws = FALSE; use coef(), ",
         "fitted() or diff_summary() instead.", call. = FALSE)
  p <- fit$p; K <- fit$K
  nsave <- dim(fit[[use]])[4]
  if (use == "C_save") {
    Omega <- fit$C_save
  } else {
    # covariance -> precision per draw per group
    Omega <- array(0, dim = c(p, p, K, nsave))
    for (s in seq_len(nsave)) {
      for (k in seq_len(K)) {
        Omega[, , k, s] <- solve(fit$Sig_save[, , k, s])
      }
    }
  }
  pcor <- array(0, dim = c(p, p, K, nsave))
  for (s in seq_len(nsave)) {
    for (k in seq_len(K)) {
      Om <- Omega[, , k, s]
      d <- sqrt(pmax(diag(Om), .Machine$double.eps))
      R <- -Om / outer(d, d)
      diag(R) <- 1
      # force symmetry for numerical noise
      pcor[, , k, s] <- 0.5 * (R + t(R))
    }
  }
  pcor
}

#' Posterior credible intervals for partial correlations
#'
#' Computes pointwise posterior credible intervals for each partial
#' correlation entry across groups.
#'
#' @param pcor A numeric array of dimension \code{[p, p, K, nsave]},
#'   as returned by \code{\link{posterior_pcor}}.
#' @param probs Numeric vector of length 2 giving the lower and upper
#'   quantile probabilities. Default \code{c(0.025, 0.975)} for 95\%
#'   credible intervals.
#'
#' @return A list with four components:
#'   \describe{
#'     \item{\code{lower}}{Numeric array \code{[p, p, K]}; lower quantile
#'       of each partial correlation.}
#'     \item{\code{median}}{Numeric array \code{[p, p, K]}; posterior median
#'       of each partial correlation.}
#'     \item{\code{upper}}{Numeric array \code{[p, p, K]}; upper quantile
#'       of each partial correlation.}
#'     \item{\code{probs}}{The quantile probabilities used.}
#'   }
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' pcor_draws <- posterior_pcor(fit)
#' ci <- posterior_ci(pcor_draws)
#' # 95% CI for edge (1,2) in group 1
#' cat(sprintf("pcor[1,2]: %.3f (%.3f, %.3f)\n",
#'     ci$median[1,2,1], ci$lower[1,2,1], ci$upper[1,2,1]))
#'
#' @seealso [posterior_pcor()], [diff_prob_pcor()]
#' @export
posterior_ci <- function(pcor, probs = c(0.025, 0.975)) {
  stopifnot(length(dim(pcor)) == 4)
  p <- dim(pcor)[1]; K <- dim(pcor)[3]; nsave <- dim(pcor)[4]

  qlo <- array(NA_real_, dim = c(p, p, K))
  q50 <- array(NA_real_, dim = c(p, p, K))
  qhi <- array(NA_real_, dim = c(p, p, K))

  for (k in seq_len(K)) {
    for (i in seq_len(p)) for (j in seq_len(p)) {
      x <- pcor[i, j, k, ]
      qs <- stats::quantile(x, probs = c(probs[1], 0.5, probs[2]), names = FALSE, na.rm = TRUE)
      qlo[i, j, k] <- qs[1]
      q50[i, j, k] <- qs[2]
      qhi[i, j, k] <- qs[3]
    }
  }
  list(lower = qlo, median = q50, upper = qhi, probs = probs)
}

#' Posterior probability of differential partial correlation (K=2)
#'
#' For two-group models, computes the posterior probability that the partial
#' correlation between each pair of variables differs between groups by
#' more than a threshold \code{delta}.
#'
#' @param pcor A numeric array of dimension \code{[p, p, 2, nsave]},
#'   as returned by \code{\link{posterior_pcor}} with \code{K = 2}.
#' @param delta Numeric; practical difference threshold. Default 0 (any
#'   difference). Set to e.g. 0.1 to detect meaningful differences only.
#'
#' @return A symmetric numeric matrix \code{[p, p]} where entry \code{[i,j]}
#'   is the posterior probability that
#'   \eqn{|\rho_{ij}^{(1)} - \rho_{ij}^{(2)}| > \delta}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' pcor_draws <- posterior_pcor(fit)
#' diff_prob <- diff_prob_pcor(pcor_draws, delta = 0.1)
#' # Edges with high probability of differential partial correlation
#' sum(diff_prob[upper.tri(diff_prob)] > 0.5)
#'
#' @seealso [posterior_pcor()], [pip_diff_edge()]
#' @export
diff_prob_pcor <- function(pcor, delta = 0) {
  stopifnot(dim(pcor)[3] == 2L)
  d <- abs(pcor[, , 1, ] - pcor[, , 2, ])
  apply(d, c(1,2), function(x) mean(x > delta, na.rm = TRUE))
}

#' Posterior inclusion probabilities (PIP) per group
#'
#' Computes the posterior inclusion probability for each edge in each group.
#' The PIP is the fraction of saved MCMC iterations in which each edge was
#' included in the graph.
#'
#' @param fit A \code{multiggm_fit} object (single chain),
#'   \code{multiggm_fit_list} object, or \code{multiggm_platform_fit} object.
#' @param chain \code{NULL} (default) pools the draws of all chains of a
#'   \code{multiggm_fit_list}; an integer vector selects chains.
#' @param platform Integer; which platform to extract (required for
#'   \code{multiggm_platform_fit} objects). Ignored for single-platform fits.
#'
#' @return A numeric array of dimension \code{[p, p, K]}. Entry
#'   \code{pip[i, j, k]} is the posterior probability that edge (i,j) is
#'   present in the graph of group \code{k}. Values range from 0 to 1.
#'   The matrix is symmetric for each group (\code{pip[i,j,k] = pip[j,i,k]}).
#'   Diagonal entries are 0.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' pip <- pip_edges(fit)
#' # Edges with PIP >= 0.5 in group 1
#' sum(pip[,,1][upper.tri(pip[,,1])] >= 0.5)
#'
#' @seealso [pip_diff_edge()], [confusion_at_threshold()], [roc_auc()]
#' @export
#'
pip_edges <- function(fit, chain = NULL, platform = NULL) {
  if (inherits(fit, "multiggm_platform_fit")) {
    if (is.null(platform)) stop("For multiggm_platform_fit, specify 'platform' argument.")
    plat <- fit$platforms[[platform]]
    A <- plat$adj_save
    return(apply(A, c(1,2,3), mean, na.rm = TRUE))
  }
  chains <- .fit_chains(fit, chain)
  if (!.all_have_draws(chains, "pip_edges()")) {
    acc <- .pool_acc(lapply(chains, .chain_edge_acc))
    return(.acc_to_array(acc$group$n_struct / acc$n, acc$p))
  }
  .pooled_mean(chains, "pip_edges()", function(A) A)
}

# TRUE if every chain saved edge-indicator draws; FALSE if some only have the
# running sums of save_draws = FALSE. Errors for fits without indicators.
.all_have_draws <- function(chains, fn) {
  for (ch in chains) .require_indicators(ch, fn)
  all(vapply(chains, function(ch) !is.null(ch$adj_save), logical(1)))
}

# Columns of a [p(p-1)/2, K] matrix (upper.tri order) as a symmetric [p, p, K] array.
.acc_to_array <- function(v, p) {
  K <- ncol(v)
  out <- array(0, c(p, p, K))
  for (k in seq_len(K)) {
    m <- matrix(0, p, p); m[upper.tri(m)] <- v[, k]
    out[, , k] <- m + t(m)
  }
  out
}

# Average of f(adj_save) over all draws of the given chains.
.pooled_mean <- function(chains, fn, f) {
  tot <- NULL; n <- 0
  for (ch in chains) {
    .require_indicators(ch, fn)
    x <- f(ch$adj_save)
    s <- apply(x, seq_len(length(dim(x)) - 1L), sum)
    tot <- if (is.null(tot)) s else tot + s
    n <- n + dim(x)[length(dim(x))]
  }
  tot / n
}

#' Posterior probability of differential edge inclusion (K=2)
#'
#' For two-group models, computes the posterior probability that each edge
#' has different inclusion status between the two groups (i.e., present in
#' one group but absent in the other).
#'
#' @param fit A \code{multiggm_fit} object (single chain) or
#'   \code{multiggm_fit_list} object. Must have \code{K = 2}.
#' @param chain \code{NULL} (default) pools the draws of all chains of a
#'   \code{multiggm_fit_list}; an integer vector selects chains.
#'
#' @return A symmetric numeric matrix \code{[p, p]} where entry \code{[i,j]}
#'   is the posterior probability that edge (i,j) has different inclusion
#'   status between groups 1 and 2. Values range from 0 to 1.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' pip_diff <- pip_diff_edge(fit)
#' sum(pip_diff[upper.tri(pip_diff)] > 0.5)
#'
#' @seealso [pip_edges()], [diff_prob_pcor()]
#' @export
pip_diff_edge <- function(fit, chain = NULL) {
  chains <- .fit_chains(fit, chain)
  if (chains[[1]]$K != 2L) stop("pip_diff_edge() needs K = 2 groups.")
  if (!.all_have_draws(chains, "pip_diff_edge()")) {
    acc <- .pool_acc(lapply(chains, .chain_edge_acc))
    return(.acc_to_array(acc$pair$n_struct / acc$n, acc$p)[, , 1])
  }
  .pooled_mean(chains, "pip_diff_edge()",
               function(A) array(A[, , 1, ] != A[, , 2, ], dim(A)[c(1, 2, 4)]))
}

#' Extract posterior summary of precision matrices
#'
#' Computes a summary statistic (default: mean) of the precision matrix
#' draws for each group. More flexible than \code{\link{coef.multiggm_fit}}
#' because you can specify an arbitrary summary function.
#'
#' @param fit A \code{multiggm_fit} or \code{multiggm_fit_list} object.
#' @param summary_fun Function to apply across posterior draws for each
#'   matrix entry. Default \code{mean}. Other useful choices include
#'   \code{median} or \code{sd}.
#' @param chain \code{NULL} (default) pools the draws of all chains of a
#'   \code{multiggm_fit_list}; an integer vector selects chains.
#'
#' @return A named list of K numeric matrices (each p x p). List names are
#'   \code{"Group_1"}, \code{"Group_2"}, etc.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' # Posterior mean
#' omega_mean <- posterior_precision(fit)
#' # Posterior standard deviation
#' omega_sd <- posterior_precision(fit, summary_fun = sd)
#'
#' @seealso [coef.multiggm_fit()], [posterior_pcor()]
#' @export
posterior_precision <- function(fit, summary_fun = mean, chain = NULL) {
  chains <- .fit_chains(fit, chain)
  K <- chains[[1]]$K; p <- chains[[1]]$p
  if (any(vapply(chains, function(ch) is.null(ch$C_save), logical(1)))) {
    if (!identical(summary_fun, mean) ||
        any(vapply(chains, function(ch) is.null(ch$C_mean), logical(1))))
      stop("This fit did not save posterior draws (save_draws = FALSE); ",
           "only the posterior mean is available.", call. = FALSE)
    wts <- vapply(chains, function(ch) as.numeric(ch$n_post), numeric(1))
    C_mean <- Reduce(`+`, Map(function(ch, w) ch$C_mean * w, chains, wts)) / sum(wts)
    out <- lapply(seq_len(K), function(k) C_mean[, , k])
  } else {
    C <- if (length(chains) == 1L) chains[[1]]$C_save else
      array(unlist(lapply(chains, `[[`, "C_save"), use.names = FALSE),
            dim = c(p, p, K, sum(vapply(chains, function(ch) dim(ch$C_save)[4], numeric(1)))))
    out <- lapply(seq_len(K), function(k)
      apply(C[, , k, , drop = FALSE], c(1, 2), summary_fun))
  }
  names(out) <- paste0("Group_", seq_len(K))
  out
}

#' Per-iteration edge counts for convergence checking
#'
#' Computes the number of edges in each group at each saved iteration.
#' Useful for trace plots of model complexity to assess MCMC convergence.
#'
#' @param fit A \code{multiggm_fit} or \code{multiggm_fit_list} object.
#' @param chain Integer; which chain to use if \code{fit} is a
#'   \code{multiggm_fit_list}. Default 1. Unlike the posterior summaries,
#'   edge counts are a per-chain trace, so chains are not pooled; use
#'   \code{plot(fit_list, type = "trace_edges")} to overlay chains.
#'
#' @return An integer matrix with \code{nsave} rows and \code{K} columns.
#'   Entry \code{[s, k]} is the number of edges (upper triangle only) in
#'   group \code{k} at saved iteration \code{s}. Column names are
#'   \code{"Group_1"}, \code{"Group_2"}, etc.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' ec <- edge_counts(fit)
#' plot(ec[,1], type = "l", ylab = "Edges", main = "Group 1 edge count trace")
#'
#' @seealso [plot_trace()], [pip_edges()]
#' @export
edge_counts <- function(fit, chain = 1L) {
  if (inherits(fit, "multiggm_fit_list")) fit <- fit$chains[[chain]]
  .require_indicators(fit, "edge_counts()")
  if (is.null(fit$adj_save))
    stop("edge_counts() needs saved draws; refit with save_draws = TRUE.", call. = FALSE)
  A <- fit$adj_save  # [p, p, K, nsave]
  nsave <- dim(A)[4]
  K <- dim(A)[3]
  p <- dim(A)[1]

  ec <- matrix(0L, nsave, K)
  for (s in seq_len(nsave)) {
    for (k in seq_len(K)) {
      # Count upper triangle edges
      ec[s, k] <- sum(A[, , k, s][upper.tri(A[, , k, s])] > 0)
    }
  }
  colnames(ec) <- paste0("Group_", seq_len(K))
  ec
}

