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

#' Extract partial correlations from posterior precision samples
#'
#' Convenience helper to convert posterior samples of precision matrices into
#' partial correlations. Supports either:
#' * a 3D array of dimension p x p x M, or
#' * a list of p x p matrices.
#'
#' @param Omega_samples Precision samples.
#' @return Partial-correlation samples in the same container type as input.
#' @export
posterior_pcor <- function(Omega_samples) {
  if (is.list(Omega_samples)) {
    return(lapply(Omega_samples, precision_to_pcor))
  }
  if (is.array(Omega_samples) && length(dim(Omega_samples)) == 3) {
    p <- dim(Omega_samples)[1]
    M <- dim(Omega_samples)[3]
    out <- array(NA_real_, dim = c(p, p, M))
    for (m in seq_len(M)) out[, , m] <- precision_to_pcor(Omega_samples[, , m])
    return(out)
  }
  stop("Omega_samples must be a list of matrices or a 3D array (p x p x M).")
}
