# Functions shared with other multiple-GGM packages (e.g. spikeyglass) live in
# ggmtools. They are re-exported so that users can call them from multiGGMr,
# and so that attaching both packages gives the same objects rather than
# conflicting versions.

#' @importFrom ggmtools precision_to_pcor
#' @export
ggmtools::precision_to_pcor

#' @importFrom ggmtools confusion_at_threshold
#' @export
ggmtools::confusion_at_threshold

#' @importFrom ggmtools roc_auc
#' @export
ggmtools::roc_auc

#' @importFrom ggmtools plot_roc
#' @export
ggmtools::plot_roc

#' @importFrom ggmtools pip_table
#' @export
ggmtools::pip_table

#' @importFrom ggmtools pcor_estimate
#' @export
ggmtools::pcor_estimate

#' @importFrom ggmtools edge_index
#' @export
ggmtools::edge_index
