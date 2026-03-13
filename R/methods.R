#' @export
print.multiggm_fit <- function(x, ...) {
  cat("<multiggm_fit>\n")
  cat("  K groups:", x$K, "\n")
  cat("  p nodes :", x$p, "\n")
  nsave <- dim(x$adj_save)[4]
  cat("  Posterior draws:", nsave, "\n")
  invisible(x)
}

#' @export
print.multiggm_fit_list <- function(x, ...) {
  cat("<multiggm_fit_list>\n")
  cat("  Chains:", x$nchains, "\n")
  cat("  K groups:", x$K, "\n")
  invisible(x)
}

#' Summarize a multiggm_fit object
#'
#' Displays MCMC diagnostics, acceptance rates, edge counts, and graph
#' similarity (theta) estimates.
#'
#' @param object A \code{multiggm_fit} object.
#' @param pip_threshold Threshold for counting selected edges. Default 0.5.
#' @param ... Ignored.
#' @return A \code{summary.multiggm_fit} object (printed invisibly).
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
#' Returns a list of K posterior mean precision matrices.
#'
#' @param object A \code{multiggm_fit} object.
#' @param ... Ignored.
#' @return A list of K numeric matrices (p x p).
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
#' Returns a list of K posterior mean partial correlation matrices.
#'
#' @param object A \code{multiggm_fit} object.
#' @param ... Ignored.
#' @return A list of K numeric matrices (p x p).
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

`%||%` <- function(a, b) if (!is.null(a)) a else b
