#' Shared-format accessors for multiGGMr fits
#'
#' Methods for the \pkg{ggmtools} generics [ggmtools::pip_table()] and
#' [ggmtools::pcor_estimate()], so that downstream code written for those
#' generics works on multiGGMr fits as well as on fits from other packages
#' (e.g. \pkg{spikeyglass}).
#'
#' \code{pip_table()} returns each group's own posterior edge-inclusion
#' probabilities (from [pip_edges()]) in long format. It needs edge
#' indicators, so it stops for mGHS fits. \code{pcor_estimate()} returns the
#' posterior mean partial correlations (as [fitted()]). For multi-chain fits
#' all chains are pooled unless \code{chain} selects one, as for
#' [pip_edges()]; for multi-platform fits, choose \code{platform}.
#'
#' multiGGMr fits do not store variable names, so the nodes are labelled
#' \code{V1, ..., Vp} unless \code{varnames} is given.
#'
#' @param fit A \code{multiggm_fit}, \code{multiggm_fit_list} or
#'   \code{multiggm_platform_fit}.
#' @param chain For a \code{multiggm_fit_list}: \code{NULL} (default) pools
#'   all chains; an integer selects chains.
#' @param platform Platform to use for a \code{multiggm_platform_fit}.
#' @param varnames Optional variable names.
#' @param ... Ignored.
#'
#' @return \code{pip_table()}: a data frame with columns \code{i}, \code{j},
#'   \code{node1}, \code{node2}, \code{group}, \code{pip}.
#'   \code{pcor_estimate()}: a list of K partial-correlation matrices named
#'   \code{Group_1}, ..., \code{Group_K}.
#'
#' @name multiggm-accessors
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 100, nsave = 50)
#' head(pip_table(fit))
#' pcor_estimate(fit)$Group_1[1:3, 1:3]
NULL

#' @rdname multiggm-accessors
#' @export
pip_table.multiggm_fit <- function(fit, varnames = NULL, ...) {
  .pip_array_to_table(pip_edges(fit), varnames)
}

#' @rdname multiggm-accessors
#' @export
pip_table.multiggm_fit_list <- function(fit, chain = NULL, varnames = NULL, ...) {
  .pip_array_to_table(pip_edges(fit, chain = chain), varnames)
}

#' @rdname multiggm-accessors
#' @export
pip_table.multiggm_platform_fit <- function(fit, platform, varnames = NULL, ...) {
  .pip_array_to_table(pip_edges(fit, platform = platform), varnames)
}

#' @rdname multiggm-accessors
#' @export
pcor_estimate.multiggm_fit <- function(fit, varnames = NULL, ...) {
  .name_pcor(fitted(fit), varnames)
}

#' @rdname multiggm-accessors
#' @export
pcor_estimate.multiggm_fit_list <- function(fit, chain = NULL, varnames = NULL, ...) {
  if (!is.null(chain)) {
    fit$chains <- fit$chains[chain]
    if (length(fit$chains) == 1L) return(.name_pcor(fitted(fit$chains[[1]]), varnames))
  }
  .name_pcor(fitted(fit), varnames)
}

#' @rdname multiggm-accessors
#' @export
pcor_estimate.multiggm_platform_fit <- function(fit, platform, varnames = NULL, ...) {
  .name_pcor(fitted(fit, platform = platform), varnames)
}

.pip_array_to_table <- function(pip, varnames) {
  K <- dim(pip)[3]
  ggmtools::as_pip_table(lapply(seq_len(K), function(k) pip[, , k]),
                         varnames = varnames,
                         groupnames = paste0("Group_", seq_len(K)))
}

.name_pcor <- function(pc, varnames) {
  p <- nrow(pc[[1]])
  if (is.null(varnames)) varnames <- paste0("V", seq_len(p))
  lapply(pc, function(m) { dimnames(m) <- list(varnames, varnames); m })
}
