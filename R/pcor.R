#' Convert precision matrix to partial correlation matrix
#'
#' Given a precision matrix \eqn{\Omega}, returns the partial correlation matrix
#' \eqn{P} with entries \eqn{P_{ij} = -\Omega_{ij}/\sqrt{\Omega_{ii}\,\Omega_{jj}}}
#' and ones on the diagonal. Partial correlations measure the linear association
#' between variables \code{i} and \code{j} after conditioning on all other
#' variables.
#'
#' @param Omega A square numeric matrix (precision / inverse covariance).
#'   Must have positive diagonal entries.
#'
#' @return A symmetric numeric matrix of the same dimension as \code{Omega},
#'   with diagonal entries equal to 1 and off-diagonal entries in \code{[-1, 1]}.
#'
#' @examples
#' Omega <- matrix(c(2, -0.5, -0.5, 1.5), 2, 2)
#' precision_to_pcor(Omega)
#'
#' @seealso [posterior_pcor()], [fitted.multiggm_fit()]
#' @export
precision_to_pcor <- function(Omega) {
  Omega <- as.matrix(Omega)
  if (nrow(Omega) != ncol(Omega)) stop("Omega must be square.")
  d <- diag(Omega)
  if (any(d <= 0)) stop("Omega must have positive diagonal to form partial correlations.")
  denom <- sqrt(outer(d, d, `*`))
  P <- -Omega / denom
  diag(P) <- 1
  P
}

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
#' @param chain Integer; which chain to use if \code{fit} is a
#'   \code{multiggm_fit_list}. Default 1.
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
posterior_pcor <- function(fit, use = c("C_save", "Sig_save"), chain = 1L) {
  use <- match.arg(use)
  if (inherits(fit, "multiggm_fit_list")) fit <- fit$chains[[chain]]
  if (!inherits(fit, "multiggm_fit")) stop("fit must be multiggm_fit or multiggm_fit_list")

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
      denom <- outer(d, d)
      R <- -Om / denom
      diag(R) <- 1
      # force symmetry for numerical noise
      R <- 0.5 * (R + t(R))
      pcor[, , k, s] <- R
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
#' @param chain Integer; which chain to use if \code{fit} is a
#'   \code{multiggm_fit_list}. Default 1.
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
pip_edges <- function(fit, chain = 1L, platform = NULL) {
  if (inherits(fit, "multiggm_platform_fit")) {
    if (is.null(platform)) stop("For multiggm_platform_fit, specify 'platform' argument.")
    plat <- fit$platforms[[platform]]
    A <- plat$adj_save
    return(apply(A, c(1,2,3), mean, na.rm = TRUE))
  }
  if (inherits(fit, "multiggm_fit_list")) fit <- fit$chains[[chain]]
  A <- fit$adj_save
  # A is [p,p,K,nsave] integer
  apply(A, c(1,2,3), mean, na.rm = TRUE)
}

#' Posterior probability of differential edge inclusion (K=2)
#'
#' For two-group models, computes the posterior probability that each edge
#' has different inclusion status between the two groups (i.e., present in
#' one group but absent in the other).
#'
#' @param fit A \code{multiggm_fit} object (single chain) or
#'   \code{multiggm_fit_list} object. Must have \code{K = 2}.
#' @param chain Integer; which chain to use if \code{fit} is a
#'   \code{multiggm_fit_list}. Default 1.
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
pip_diff_edge <- function(fit, chain = 1L) {
  if (inherits(fit, "multiggm_fit_list")) fit <- fit$chains[[chain]]
  A <- fit$adj_save
  stopifnot(dim(A)[3] == 2L)
  d <- (A[, , 1, ] != A[, , 2, ])
  apply(d, c(1,2), mean, na.rm = TRUE)
}

#' Confusion metrics at a PIP threshold
#'
#' Computes true/false positive/negative counts and rates by comparing a
#' continuous score matrix (e.g., PIP) against a binary ground truth
#' adjacency matrix at a given threshold. Only the upper triangle is used.
#'
#' @param score_mat Numeric matrix \code{[p, p]}; continuous edge scores
#'   (e.g., posterior inclusion probabilities from \code{\link{pip_edges}}).
#' @param truth_mat Numeric matrix \code{[p, p]}; binary ground truth
#'   adjacency (1 = edge present, 0 = absent).
#' @param thr Numeric; threshold for edge selection. An edge is selected if
#'   \code{score_mat[i,j] >= thr}.
#' @param upper_only Logical; if \code{TRUE} (default), only the upper
#'   triangle is evaluated (appropriate for undirected graphs).
#'
#' @return A named numeric vector with components:
#'   \describe{
#'     \item{\code{TP}}{True positives (correctly selected edges).}
#'     \item{\code{FP}}{False positives (incorrectly selected non-edges).}
#'     \item{\code{TN}}{True negatives (correctly excluded non-edges).}
#'     \item{\code{FN}}{False negatives (missed true edges).}
#'     \item{\code{TPR}}{True positive rate (sensitivity / recall).}
#'     \item{\code{FPR}}{False positive rate (1 - specificity).}
#'   }
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' pip <- pip_edges(fit)
#' confusion_at_threshold(pip[,,1], sim$adj_list[[1]], thr = 0.5)
#'
#' @seealso [roc_auc()], [pip_edges()]
#' @export
confusion_at_threshold <- function(score_mat, truth_mat, thr, upper_only = TRUE) {
  stopifnot(all(dim(score_mat) == dim(truth_mat)))
  p <- nrow(score_mat)

  idx <- which(upper.tri(score_mat), arr.ind = TRUE)
  s <- score_mat[idx]
  y <- truth_mat[idx] > 0

  pred <- s >= thr
  TP <- sum(pred & y)
  FP <- sum(pred & !y)
  TN <- sum(!pred & !y)
  FN <- sum(!pred & y)

  TPR <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  FPR <- if ((FP + TN) > 0) FP / (FP + TN) else NA_real_
  c(TP=TP, FP=FP, TN=TN, FN=FN, TPR=TPR, FPR=FPR)
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
#' @param chain Integer; which chain to use if \code{fit} is a
#'   \code{multiggm_fit_list}. Default 1.
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
posterior_precision <- function(fit, summary_fun = mean, chain = 1L) {
  if (inherits(fit, "multiggm_fit_list")) fit <- fit$chains[[chain]]
  K <- fit$K
  out <- vector("list", K)
  for (k in seq_len(K)) {
    out[[k]] <- apply(fit$C_save[, , k, , drop = FALSE], c(1, 2), summary_fun)
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
#'   \code{multiggm_fit_list}. Default 1.
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

#' ROC curve and AUC from continuous scores
#'
#' Computes the receiver operating characteristic (ROC) curve and area
#' under the curve (AUC) by comparing a continuous score matrix (e.g., PIP)
#' against a binary ground truth adjacency matrix. Only the upper triangle
#' is used.
#'
#' @param score_mat Numeric matrix \code{[p, p]}; continuous edge scores
#'   (e.g., posterior inclusion probabilities from \code{\link{pip_edges}}).
#' @param truth_mat Numeric matrix \code{[p, p]}; binary ground truth
#'   adjacency (1 = edge present, 0 = absent).
#'
#' @return A list with three components:
#'   \describe{
#'     \item{\code{FPR}}{Numeric vector; false positive rates at each
#'       threshold (starts at 0, ends at 1).}
#'     \item{\code{TPR}}{Numeric vector; true positive rates at each
#'       threshold (starts at 0, ends at 1).}
#'     \item{\code{auc}}{Numeric scalar; area under the ROC curve.
#'       Values near 1 indicate excellent discrimination.}
#'   }
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' pip <- pip_edges(fit)
#' roc <- roc_auc(pip[,,1], sim$adj_list[[1]])
#' cat("AUC:", roc$auc, "\n")
#'
#' @seealso [plot_roc()], [confusion_at_threshold()], [pip_edges()]
#' @export
roc_auc <- function(score_mat, truth_mat) {
  stopifnot(all(dim(score_mat) == dim(truth_mat)))
  idx <- which(upper.tri(score_mat), arr.ind = TRUE)
  s <- score_mat[idx]
  y <- truth_mat[idx] > 0

  ord <- order(s, decreasing = TRUE)
  s <- s[ord]; y <- y[ord]

  P <- sum(y); N <- sum(!y)
  if (P == 0 || N == 0) stop("Need both positives and negatives in truth_mat.")

  TP <- cumsum(y)
  FP <- cumsum(!y)
  TPR <- TP / P
  FPR <- FP / N

  # trapezoidal AUC in FPR-TPR space
  auc <- sum(diff(c(0, FPR)) * (head(c(0, TPR), -1) + tail(c(0, TPR), -1)) / 2)

  list(FPR = c(0, FPR), TPR = c(0, TPR), auc = auc)
}
