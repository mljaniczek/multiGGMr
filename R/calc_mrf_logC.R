#' Log normalizing constant for the MRF prior (per edge)
#'
#' Computes log C(Theta, nu) where
#' C(Theta, nu) = sum_{g in {0,1}^K} exp( nu * sum(g) + g' Theta g ).
#'
#' This is the computational hot spot in the reference MATLAB implementation
#' and is implemented in C++ via RcppArmadillo.
#'
#' @param Theta K x K numeric matrix (graph similarity parameters).
#' @param nu Numeric vector of length 1 or more (edge-specific log-odds).
#' @return Numeric vector logC of the same length as nu.
#' @examples
#' Theta <- matrix(0, 3, 3)
#' nu <- c(-1, 0, 1)
#' calc_mrf_logC(Theta, nu)
#' @export
calc_mrf_logC <- function(Theta, nu) {
  Theta <- as.matrix(Theta)
  nu <- as.numeric(nu)
  calc_mrf_logC_cpp(Theta, nu)
}
