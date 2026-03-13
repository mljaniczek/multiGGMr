#' Rescale a symmetric matrix toward positive definiteness
#'
#' Port of the helper in the MATLAB reference code. For each row, rescales
#' off-diagonal entries by \code{(denom_factor * row L1 norm of off-diagonals)},
#' then enforces unit diagonal and symmetry. This makes the matrix
#' diagonally dominant (hence positive definite) when \code{denom_factor > 1}.
#'
#' @param A Numeric square matrix with unit diagonal.
#' @param denom_factor Positive scalar controlling shrinkage of off-diagonal
#'   entries. Values > 1 shrink off-diagonals more aggressively. Default 1.
#'
#' @return A symmetric numeric matrix with unit diagonal, where each row's
#'   off-diagonal entries have been rescaled so their absolute sum is at
#'   most \code{1 / denom_factor}.
#'
#' @examples
#' A <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
#' fix_matrix(A, denom_factor = 1.5)
#'
#' @export
fix_matrix <- function(A, denom_factor = 1) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  stopifnot(is.numeric(denom_factor), length(denom_factor) == 1, denom_factor > 0)

  p <- nrow(A)
  A <- as.matrix(A)

  for (cur_row in seq_len(p)) {
    cur_sum <- sum(abs(A[cur_row, ])) - 1
    if (!isTRUE(all.equal(cur_sum, 0))) {
      A[cur_row, ] <- A[cur_row, ] / (denom_factor * cur_sum)
    }
    A[cur_row, cur_row] <- 1
  }

  # symmetrize
  (A + t(A)) / 2
}
