#' Ensure a symmetric positive-definite-ish correlation/precision template
#'
#' Port of the helper in the MATLAB reference code. For each row, rescales
#' off-diagonal entries by the (denom_factor * row L1 norm of off-diagonals),
#' then enforces unit diagonal and symmetry.
#'
#' @param A Numeric square matrix.
#' @param denom_factor Positive scalar. Values > 1 shrink off-diagonals more.
#' @return A symmetric matrix with unit diagonal.
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
