#' Convert precision matrix to partial correlation matrix
#'
#' Given a precision matrix \eqn{\Omega}, returns the partial correlation matrix
#' \eqn{P} with entries \eqn{p_{ij} = -\omega_{ij}/\sqrt{\omega_{ii}\,\omega_{jj}}}
#' and ones on the diagonal.
#'
#' @param Omega A square numeric matrix (precision).
#' @return A numeric matrix of partial correlations with diagonal 1.
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
#' @param fit multiggm_fit (single chain) or multiggm_fit_list
#' @param use which posterior array to use: "C_save" (precision) or "Sig_save" (covariance -> invert)
#' @param chain which chain (if fit_list)
#' @note for each draw and group k, partial correlation is: rho_ij^{(k)} = -Omega_{ij}^{(k)}/sqrt{Omega_{ii}^{(k)}Omega_{jj}^{(k)}}
#' @return array [p, p, K, nsave] of partial correlation draws
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
#' @param pcor array [p,p,K,nsave]
#' @param probs numeric length-2, e.g. c(0.025, 0.975)
#' @return list with lower/median/upper arrays [p,p,K]
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
#' @param pcor array [p,p,2,nsave]
#' @param delta practical difference threshold (0 = any difference)
#' @return matrix [p,p] of posterior probabilities
#' @export
diff_prob_pcor <- function(pcor, delta = 0) {
  stopifnot(dim(pcor)[3] == 2L)
  d <- abs(pcor[, , 1, ] - pcor[, , 2, ])
  apply(d, c(1,2), function(x) mean(x > delta, na.rm = TRUE))
}

#' Posterior inclusion probabilities (PIP) per group
#' @param fit a multi_ggm_fit_list object
#' @return array [p,p,K]
#' @export
pip_edges <- function(fit, chain = 1L) {
  if (inherits(fit, "multiggm_fit_list")) fit <- fit$chains[[chain]]
  A <- fit$adj_save
  # A is [p,p,K,nsave] integer
  apply(A, c(1,2,3), mean, na.rm = TRUE)
}

#' Posterior inclusion probabilities (PIP) per group differential edges
#' @return array [p,p,K]
#' @export
pip_diff_edge <- function(fit, chain = 1L) {
  if (inherits(fit, "multiggm_fit_list")) fit <- fit$chains[[chain]]
  A <- fit$adj_save
  stopifnot(dim(A)[3] == 2L)
  d <- (A[, , 1, ] != A[, , 2, ])
  apply(d, c(1,2), mean, na.rm = TRUE)
}

#' Confusion metrics at a threshold
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
#' Returns a list of K matrices summarizing the posterior precision draws.
#'
#' @param fit A \code{multiggm_fit} or \code{multiggm_fit_list} object.
#' @param summary_fun Function to apply across posterior draws. Default \code{mean}.
#' @param chain Which chain (if \code{fit_list}). Default 1.
#' @return A list of K numeric matrices (p x p).
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
#' Useful for trace plots of model complexity.
#'
#' @param fit A \code{multiggm_fit} or \code{multiggm_fit_list} object.
#' @param chain Which chain (if \code{fit_list}). Default 1.
#' @return A matrix with \code{nsave} rows and K columns.
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

#' ROC and AUC from continuous scores
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




