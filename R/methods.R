#' Print a multiggm_fit object
#'
#' Displays a concise summary of a fitted multi-GGM model.
#'
#' @param x A \code{multiggm_fit} object returned by \code{\link{multiggm_mcmc}}.
#' @param ... Ignored.
#' @return Invisible \code{x}.
#' @export
print.multiggm_fit <- function(x, ...) {
  cat("<multiggm_fit>\n")
  cat("  K groups:", x$K, "\n")
  cat("  p nodes :", x$p, "\n")
  nsave <- dim(x$adj_save)[4]
  cat("  Posterior draws:", nsave, "\n")
  invisible(x)
}

#' Print a multiggm_fit_list object
#'
#' Displays a concise summary of a multi-chain multi-GGM fit.
#'
#' @param x A \code{multiggm_fit_list} object returned by
#'   \code{\link{multiggm_mcmc}} with \code{nchains > 1}.
#' @param ... Ignored.
#' @return Invisible \code{x}.
#' @export
print.multiggm_fit_list <- function(x, ...) {
  cat("<multiggm_fit_list>\n")
  cat("  Chains:", x$nchains, "\n")
  cat("  K groups:", x$K, "\n")
  invisible(x)
}

#' Summarize a multiggm_fit object
#'
#' Displays MCMC diagnostics including acceptance rates, edge counts at a
#' given PIP threshold, and graph similarity (theta) estimates.
#'
#' @param object A \code{multiggm_fit} object returned by
#'   \code{\link{multiggm_mcmc}}.
#' @param pip_threshold Numeric threshold for counting selected edges via
#'   posterior inclusion probability (PIP). Default 0.5.
#' @param ... Ignored.
#'
#' @return An object of class \code{"summary.multiggm_fit"} (printed
#'   invisibly) with components:
#'   \describe{
#'     \item{\code{K}}{Integer; number of groups.}
#'     \item{\code{p}}{Integer; number of variables.}
#'     \item{\code{nsave}}{Integer; number of saved posterior draws.}
#'     \item{\code{ar_gamma}}{Numeric; mean acceptance rate for between-model
#'       (spike-slab toggle) moves on theta.}
#'     \item{\code{ar_theta}}{Numeric; mean acceptance rate for within-model
#'       (slab-to-slab) moves on theta.}
#'     \item{\code{ar_nu}}{Numeric; mean acceptance rate for nu updates.}
#'     \item{\code{edge_counts}}{Integer vector of length K; number of selected
#'       edges per group at the given PIP threshold.}
#'     \item{\code{pip_threshold}}{The PIP threshold used.}
#'     \item{\code{theta_mean}}{Numeric matrix \code{[K, K]}; posterior mean of
#'       theta (graph similarity) for each pair of groups. Upper triangle only.}
#'     \item{\code{theta_nonzero_frac}}{Numeric matrix \code{[K, K]}; posterior
#'       probability that theta > 0 for each pair. Upper triangle only.}
#'     \item{\code{hyper}}{Named list of hyperparameters used in the fit.}
#'   }
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' summary(fit)
#' summary(fit, pip_threshold = 0.3)
#'
#' @seealso [multiggm_mcmc()], [pip_edges()]
#' @export
summary.multiggm_fit <- function(object, pip_threshold = 0.5, ...) {
  K <- object$K
  p <- object$p
  nsave <- dim(object$adj_save)[4]

  # Acceptance rates (upper triangle only for symmetric quantities)
  ar_gamma_mean <- if (K > 1) mean(object$ar_gamma[upper.tri(object$ar_gamma)]) else NA_real_
  ar_theta_mean <- if (K > 1) mean(object$ar_theta[upper.tri(object$ar_theta)]) else NA_real_
  ar_nu_mean    <- mean(object$ar_nu[upper.tri(object$ar_nu)])

  # Edge counts at PIP threshold
  pip <- pip_edges(object)
  edge_counts <- integer(K)
  for (k in seq_len(K)) {
    pip_k <- pip[, , k]
    edge_counts[k] <- sum(pip_k[upper.tri(pip_k)] >= pip_threshold)
  }

  # Theta summary
  theta_mean <- NULL
  theta_nonzero_frac <- NULL
  if (K > 1) {
    theta_mean <- matrix(0, K, K)
    theta_nonzero_frac <- matrix(0, K, K)
    for (k in 1:(K - 1)) {
      for (m in (k + 1):K) {
        vals <- object$Theta_save[k, m, ]
        theta_mean[k, m] <- mean(vals)
        theta_mean[m, k] <- theta_mean[k, m]
        theta_nonzero_frac[k, m] <- mean(vals > 0)
        theta_nonzero_frac[m, k] <- theta_nonzero_frac[k, m]
      }
    }
  }

  out <- list(
    K = K,
    p = p,
    nsave = nsave,
    ar_gamma = ar_gamma_mean,
    ar_theta = ar_theta_mean,
    ar_nu = ar_nu_mean,
    edge_counts = edge_counts,
    pip_threshold = pip_threshold,
    theta_mean = theta_mean,
    theta_nonzero_frac = theta_nonzero_frac,
    hyper = object$hyper
  )
  class(out) <- "summary.multiggm_fit"
  out
}

#' @export
print.summary.multiggm_fit <- function(x, ...) {
  cat("multiGGM MCMC Summary\n")
  cat("=====================\n")
  cat("Groups (K):", x$K, "  |  Nodes (p):", x$p, "  |  Posterior draws:", x$nsave, "\n\n")

  cat("Acceptance rates:\n")
  if (!is.na(x$ar_gamma))
    cat("  gamma (edge toggle):", sprintf("%.1f%%", x$ar_gamma * 100), "\n")
  if (!is.na(x$ar_theta))
    cat("  theta (within-model):", sprintf("%.1f%%", x$ar_theta * 100), "\n")
  cat("  nu (edge log-odds):", sprintf("%.1f%%", x$ar_nu * 100), "\n\n")

  cat("Selected edges (PIP >=", x$pip_threshold, "):\n")
  for (k in seq_along(x$edge_counts)) {
    cat("  Group", k, ":", x$edge_counts[k], "edges\n")
  }

  if (!is.null(x$theta_mean) && x$K > 1) {
    cat("\nGraph similarity (theta):\n")
    for (k in 1:(x$K - 1)) {
      for (m in (k + 1):x$K) {
        cat(sprintf("  theta[%d,%d]: mean = %.3f, P(nonzero) = %.1f%%\n",
                    k, m, x$theta_mean[k, m], x$theta_nonzero_frac[k, m] * 100))
      }
    }
  }

  cat("\n")
  invisible(x)
}

#' Extract posterior mean precision matrices
#'
#' Returns a list of K posterior mean precision matrices from a fitted
#' multi-GGM model. This is the \code{coef} S3 method for
#' \code{multiggm_fit} objects.
#'
#' @param object A \code{multiggm_fit} object returned by
#'   \code{\link{multiggm_mcmc}}.
#' @param ... Ignored.
#'
#' @return A named list of K numeric matrices (each p x p). Each matrix is
#'   the element-wise posterior mean of the precision matrix \eqn{\Omega_k}
#'   across all saved MCMC iterations. List names are \code{"Group_1"},
#'   \code{"Group_2"}, etc.
#'
#' @details
#' The precision matrix \eqn{\Omega_k} encodes the conditional independence
#' structure of group \code{k}: \eqn{\Omega_{ij}^{(k)} = 0} if and only if
#' variables \code{i} and \code{j} are conditionally independent in group
#' \code{k}. The posterior mean is computed from \code{object$C_save}.
#'
#' For partial correlations (standardized precision), use
#' \code{\link{fitted.multiggm_fit}} instead.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' omega_hat <- coef(fit)
#' round(omega_hat$Group_1[1:5, 1:5], 3)
#'
#' @seealso [fitted.multiggm_fit()], [posterior_precision()],
#'   [posterior_pcor()]
#' @export
coef.multiggm_fit <- function(object, ...) {
  K <- object$K
  p <- object$p
  out <- vector("list", K)
  for (k in seq_len(K)) {
    out[[k]] <- apply(object$C_save[, , k, , drop = FALSE], c(1, 2), mean)
  }
  names(out) <- paste0("Group_", seq_len(K))
  out
}

#' Extract posterior mean partial correlation matrices
#'
#' Returns a list of K posterior mean partial correlation matrices from a
#' fitted multi-GGM model. This is the \code{fitted} S3 method for
#' \code{multiggm_fit} objects.
#'
#' @param object A \code{multiggm_fit} object returned by
#'   \code{\link{multiggm_mcmc}}.
#' @param ... Ignored.
#'
#' @return A named list of K numeric matrices (each p x p). Each matrix is
#'   the element-wise posterior mean of the partial correlation matrix
#'   \eqn{P_k} across all saved MCMC iterations, where
#'   \eqn{P_{ij}^{(k)} = -\Omega_{ij}^{(k)} / \sqrt{\Omega_{ii}^{(k)} \Omega_{jj}^{(k)}}}.
#'   Diagonal entries are 1. List names are \code{"Group_1"},
#'   \code{"Group_2"}, etc.
#'
#' @details
#' For each saved iteration, the precision matrix \eqn{\Omega_k} is
#' converted to a partial correlation matrix, then the posterior mean is
#' taken element-wise. This is computed via \code{\link{posterior_pcor}}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' pcor_hat <- fitted(fit)
#' round(pcor_hat$Group_1[1:5, 1:5], 3)
#'
#' @seealso [coef.multiggm_fit()], [posterior_pcor()],
#'   [precision_to_pcor()]
#' @export
fitted.multiggm_fit <- function(object, ...) {
  pcor <- posterior_pcor(object)
  K <- object$K
  out <- vector("list", K)
  for (k in seq_len(K)) {
    out[[k]] <- apply(pcor[, , k, , drop = FALSE], c(1, 2), mean)
  }
  names(out) <- paste0("Group_", seq_len(K))
  out
}

# ============================================================================
# Multi-platform fit S3 methods
# ============================================================================

#' Print a multiggm_platform_fit object
#'
#' Displays a concise summary of a fitted multi-platform multi-GGM model.
#'
#' @param x A \code{multiggm_platform_fit} object returned by
#'   \code{\link{multiggm_mcmc}} with \code{method = "ssvs_platform"}.
#' @param ... Ignored.
#' @return Invisible \code{x}.
#' @export
print.multiggm_platform_fit <- function(x, ...) {
  cat("<multiggm_platform_fit>\n")
  cat("  Platforms (S):", x$S, "\n")
  cat("  Groups (K):", x$K, "\n")
  cat("  Nodes per platform:", paste(x$p_vec, collapse = ", "), "\n")
  nsave <- dim(x$Phi_save)[3]
  cat("  Posterior draws:", nsave, "\n")
  invisible(x)
}

#' Summarize a multiggm_platform_fit object
#'
#' Displays per-platform MCMC diagnostics and platform-level similarity
#' estimates.
#'
#' @param object A \code{multiggm_platform_fit} object.
#' @param pip_threshold Numeric threshold for counting selected edges.
#'   Default 0.5.
#' @param ... Ignored.
#' @return An object of class \code{"summary.multiggm_platform_fit"} (printed
#'   invisibly).
#' @export
summary.multiggm_platform_fit <- function(object, pip_threshold = 0.5, ...) {
  K <- object$K
  S <- object$S
  nsave <- dim(object$Phi_save)[3]

  # Per-platform summaries
  platform_summaries <- vector("list", S)
  for (s in seq_len(S)) {
    plat <- object$platforms[[s]]
    p_s <- plat$p
    ar_gamma_mean <- if (K > 1) mean(plat$ar_gamma[upper.tri(plat$ar_gamma)]) else NA_real_
    ar_theta_mean <- if (K > 1) mean(plat$ar_theta[upper.tri(plat$ar_theta)]) else NA_real_
    ar_nu_mean <- mean(plat$ar_nu[upper.tri(plat$ar_nu)])

    # Edge counts at PIP threshold
    pip_arr <- apply(plat$adj_save, c(1, 2, 3), mean)
    edge_counts <- integer(K)
    for (k in seq_len(K)) {
      pip_k <- pip_arr[, , k]
      edge_counts[k] <- sum(pip_k[upper.tri(pip_k)] >= pip_threshold)
    }

    # Theta summary
    theta_mean <- NULL
    theta_nonzero_frac <- NULL
    if (K > 1) {
      theta_mean <- matrix(0, K, K)
      theta_nonzero_frac <- matrix(0, K, K)
      for (k in 1:(K - 1)) {
        for (m in (k + 1):K) {
          vals <- plat$Theta_save[k, m, ]
          theta_mean[k, m] <- mean(vals)
          theta_mean[m, k] <- theta_mean[k, m]
          theta_nonzero_frac[k, m] <- mean(vals > 0)
          theta_nonzero_frac[m, k] <- theta_nonzero_frac[k, m]
        }
      }
    }

    platform_summaries[[s]] <- list(
      p = p_s, ar_gamma = ar_gamma_mean, ar_theta = ar_theta_mean,
      ar_nu = ar_nu_mean, edge_counts = edge_counts,
      theta_mean = theta_mean, theta_nonzero_frac = theta_nonzero_frac
    )
  }

  # Platform-level Phi summary
  phi_mean <- NULL
  phi_nonzero_frac <- NULL
  ar_phi_between_mean <- NA_real_
  ar_phi_within_mean <- NA_real_
  ar_w_mean <- NA_real_

  if (S > 1) {
    phi_mean <- matrix(0, S, S)
    phi_nonzero_frac <- matrix(0, S, S)
    for (s1 in 1:(S - 1)) {
      for (s2 in (s1 + 1):S) {
        vals <- object$Phi_save[s1, s2, ]
        phi_mean[s1, s2] <- mean(vals)
        phi_mean[s2, s1] <- phi_mean[s1, s2]
        phi_nonzero_frac[s1, s2] <- mean(vals > 0)
        phi_nonzero_frac[s2, s1] <- phi_nonzero_frac[s1, s2]
      }
    }
    ar_phi_between_mean <- mean(object$ar_phi_between[upper.tri(object$ar_phi_between)])
    ar_phi_within_mean <- mean(object$ar_phi_within[upper.tri(object$ar_phi_within)])
    if (K > 1) ar_w_mean <- mean(object$ar_w[upper.tri(object$ar_w)])
  }

  out <- list(
    K = K, S = S, p_vec = object$p_vec, nsave = nsave,
    pip_threshold = pip_threshold,
    platform_summaries = platform_summaries,
    phi_mean = phi_mean,
    phi_nonzero_frac = phi_nonzero_frac,
    ar_phi_between = ar_phi_between_mean,
    ar_phi_within = ar_phi_within_mean,
    ar_w = ar_w_mean
  )
  class(out) <- "summary.multiggm_platform_fit"
  out
}

#' @export
print.summary.multiggm_platform_fit <- function(x, ...) {
  cat("multiGGM Multi-Platform MCMC Summary\n")
  cat("=====================================\n")
  cat("Platforms (S):", x$S, "  |  Groups (K):", x$K,
      "  |  Posterior draws:", x$nsave, "\n")
  cat("Nodes per platform:", paste(x$p_vec, collapse = ", "), "\n\n")

  for (s in seq_along(x$platform_summaries)) {
    ps <- x$platform_summaries[[s]]
    cat("--- Platform", s, "(p =", ps$p, ") ---\n")
    cat("  Acceptance rates:\n")
    if (!is.na(ps$ar_gamma))
      cat("    gamma (edge toggle):", sprintf("%.1f%%", ps$ar_gamma * 100), "\n")
    if (!is.na(ps$ar_theta))
      cat("    theta (within-model):", sprintf("%.1f%%", ps$ar_theta * 100), "\n")
    cat("    nu (edge log-odds):", sprintf("%.1f%%", ps$ar_nu * 100), "\n")

    cat("  Selected edges (PIP >=", x$pip_threshold, "):\n")
    for (k in seq_along(ps$edge_counts))
      cat("    Group", k, ":", ps$edge_counts[k], "edges\n")

    if (!is.null(ps$theta_mean) && x$K > 1) {
      cat("  Graph similarity (theta):\n")
      for (k in 1:(x$K - 1))
        for (m in (k + 1):x$K)
          cat(sprintf("    theta[%d,%d]: mean = %.3f, P(nonzero) = %.1f%%\n",
                      k, m, ps$theta_mean[k, m], ps$theta_nonzero_frac[k, m] * 100))
    }
    cat("\n")
  }

  if (x$S > 1) {
    cat("--- Platform-level ---\n")
    cat("  Acceptance rates:\n")
    cat("    Phi (between-model):", sprintf("%.1f%%", x$ar_phi_between * 100), "\n")
    cat("    Phi (within-model):", sprintf("%.1f%%", x$ar_phi_within * 100), "\n")
    if (!is.na(x$ar_w))
      cat("    w (platform log-odds):", sprintf("%.1f%%", x$ar_w * 100), "\n")

    if (!is.null(x$phi_mean)) {
      cat("  Platform similarity (Phi):\n")
      for (s1 in 1:(x$S - 1))
        for (s2 in (s1 + 1):x$S)
          cat(sprintf("    Phi[%d,%d]: mean = %.3f, P(nonzero) = %.1f%%\n",
                      s1, s2, x$phi_mean[s1, s2], x$phi_nonzero_frac[s1, s2] * 100))
    }
  }

  cat("\n")
  invisible(x)
}

#' Extract posterior mean precision matrices from a multi-platform fit
#'
#' @param object A \code{multiggm_platform_fit} object.
#' @param platform Integer; which platform to extract. If NULL (default),
#'   returns all platforms.
#' @param ... Ignored.
#' @return A list of K posterior mean precision matrices for the specified
#'   platform, or a list of S such lists.
#' @export
coef.multiggm_platform_fit <- function(object, platform = NULL, ...) {
  if (!is.null(platform)) {
    plat <- object$platforms[[platform]]
    K <- object$K
    out <- vector("list", K)
    for (k in seq_len(K))
      out[[k]] <- apply(plat$C_save[, , k, , drop = FALSE], c(1, 2), mean)
    names(out) <- paste0("Group_", seq_len(K))
    return(out)
  }

  out <- vector("list", object$S)
  for (s in seq_len(object$S)) {
    out[[s]] <- coef.multiggm_platform_fit(object, platform = s)
  }
  names(out) <- paste0("Platform_", seq_len(object$S))
  out
}

#' Extract posterior mean partial correlation matrices from a multi-platform fit
#'
#' @param object A \code{multiggm_platform_fit} object.
#' @param platform Integer; which platform. If NULL, returns all.
#' @param ... Ignored.
#' @return A list of K posterior mean partial correlation matrices, or a list
#'   of S such lists.
#' @export
fitted.multiggm_platform_fit <- function(object, platform = NULL, ...) {
  if (!is.null(platform)) {
    plat <- object$platforms[[platform]]
    K <- object$K
    nsave <- dim(plat$C_save)[4]
    p_s <- plat$p
    pcor_arr <- array(NA_real_, dim = c(p_s, p_s, K, nsave))
    for (s_idx in seq_len(nsave)) {
      for (k in seq_len(K)) {
        pcor_arr[, , k, s_idx] <- precision_to_pcor(plat$C_save[, , k, s_idx])
      }
    }
    out <- vector("list", K)
    for (k in seq_len(K))
      out[[k]] <- apply(pcor_arr[, , k, , drop = FALSE], c(1, 2), mean)
    names(out) <- paste0("Group_", seq_len(K))
    return(out)
  }

  out <- vector("list", object$S)
  for (s in seq_len(object$S)) {
    out[[s]] <- fitted.multiggm_platform_fit(object, platform = s)
  }
  names(out) <- paste0("Platform_", seq_len(object$S))
  out
}


`%||%` <- function(a, b) if (!is.null(a)) a else b
