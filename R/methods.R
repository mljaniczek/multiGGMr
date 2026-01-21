#' @export
print.multiggm_fit <- function(x, ...) {
  cat("<multiggm_fit>\n")
  if (!is.null(x$K)) cat("K graphs:", x$K, "\n")
  if (!is.null(x$p)) cat("p nodes :", x$p, "\n")
  invisible(x)
}

#' @export
print.multiggm_fit_list <- function(x, ...) {
  cat("<multiggm_fit_list>\n")
  cat("Chains:", x$nchains, "\n")
  cat("K graphs:", x$K, "\n")
  invisible(x)
}

#' @export
summary.multiggm_fit <- function(object, ...) {
  # Placeholder summary, to be expanded once samples are stored.
  out <- list(call = object$call %||% NULL,
              K = object$K %||% NULL,
              p = object$p %||% NULL)
  class(out) <- "summary.multiggm_fit"
  out
}

#' @export
print.summary.multiggm_fit <- function(x, ...) {
  cat("<summary.multiggm_fit>\n")
  if (!is.null(x$K)) cat("K graphs:", x$K, "\n")
  if (!is.null(x$p)) cat("p nodes :", x$p, "\n")
  invisible(x)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
