# Posterior predictive checks: replicate data from posterior draws of the
# precision matrices and compare discrepancy statistics with the observed
# data. Runs in R after fitting, so it never touches the sampler's RNG stream.

# Sample correlations and, when n > p + 1, partial correlations of a centered
# data matrix.
.ppc_cor <- function(X, edge_stat) {
  S <- crossprod(X)
  if (edge_stat == "cor") return(stats::cov2cor(S))
  P <- -stats::cov2cor(solve(S))
  diag(P) <- 1
  P
}

# Discrepancies that depend on the data only (per group).
.ppc_group_stats <- function(R) {
  r <- R[upper.tri(R)]
  c(max_abs_cor = max(abs(r)), sum_sq_cor = sum(r^2))
}

# Realized discrepancies (Gelman, Meng & Stern 1996), which also depend on
# the precision draw: squared Mahalanobis distances y_i' Omega y_i.
.ppc_realized <- function(X, U) {
  m <- rowSums((X %*% t(U))^2)
  c(max_mahal = max(m), kurtosis = mean(m^2))
}

#' Posterior predictive checks
#'
#' Checks whether the fitted model reproduces features of the observed data.
#' For each of \code{ndraws} posterior draws of the precision matrices,
#' replicate data with the observed group sizes are drawn from the model, and
#' discrepancy statistics of the replicate data are compared with those of
#' the observed data. A posterior predictive p-value (ppp) near 0 or 1 marks
#' a feature of the data the model does not reproduce.
#'
#' @param fit A \code{multiggm_fit} or \code{multiggm_fit_list} (G-Wishart,
#'   SSVS or mGHS).
#' @param data_list The list of data matrices the model was fitted to. They
#'   are column-centered as in [multiggm_mcmc()].
#' @param ndraws Number of posterior draws to use, spread evenly over the
#'   saved draws of all chains. Default 200.
#' @param use \code{"draws"} (default) uses saved precision draws, so the
#'   fit needs \code{save_draws = TRUE} (with \code{thin} to keep it small).
#'   \code{"mean"} replicates every dataset from the posterior mean precision
#'   matrices, for fits without draws. It ignores posterior uncertainty, so
#'   its p-values are more extreme than true posterior predictive p-values;
#'   treat it as approximate.
#' @param edge_stat Edge-level statistic: \code{"pcor"} (sample partial
#'   correlation; needs \eqn{n_k > p + 1}), \code{"cor"} (sample
#'   correlation), or \code{"auto"} (default: \code{"pcor"} when every group
#'   allows it).
#' @param alpha Two-sided level used to count flagged edges in the printed
#'   summary. Default 0.05.
#' @param chain \code{NULL} (default) pools all chains.
#' @param seed Optional seed for the replicate data (sets R's global RNG).
#'
#' @details
#' Statistics, with ppp = \eqn{P(T(y^{rep}) \ge T(y) \mid y)}:
#' \describe{
#'   \item{\code{max_abs_cor}, \code{sum_sq_cor}}{Per group: the largest
#'     absolute sample correlation and the sum of squared sample correlations
#'     (overall strength of dependence). A ppp near 1 for \code{sum_sq_cor}
#'     means the model implies more dependence than observed, near 0 less
#'     (for example from over-shrinkage).}
#'   \item{\code{max_mahal}, \code{kurtosis}}{Per group, realized
#'     discrepancies computed with the same precision draw for observed and
#'     replicate data: the largest squared Mahalanobis distance
#'     \eqn{y_i' \Omega y_i} (outliers) and their mean square, Mardia's
#'     multivariate kurtosis (heavy tails). A ppp near 0 means the data have
#'     more extreme observations or heavier tails than a Gaussian model
#'     produces.}
#'   \item{\code{sum_sq_diff}}{Per pair of groups: the sum over edges of the
#'     squared difference in the edge statistic between the groups. A ppp
#'     near 0 means the observed groups differ more than the fitted model
#'     implies (differences are over-shrunk); near 1, less.}
#' }
#' Per edge, the same comparison is made for the edge statistic in each group
#' and for its difference between each pair of groups (\code{edges}).
#'
#' Posterior predictive p-values are not uniform under a correct model: they
#' concentrate around 0.5, so a correct model flags fewer edges than
#' \code{alpha} would suggest. Many flagged edges, or extreme global ppp
#' values, are evidence of misfit; few flags do not prove fit.
#'
#' @return An object of class \code{"multiggm_ppc"}: a list with
#'   \describe{
#'     \item{\code{global}}{Data frame: \code{type} (\code{"group"} or
#'       \code{"pair"}), \code{group1}, \code{group2}, \code{statistic},
#'       \code{observed} (averaged over draws for the realized statistics),
#'       \code{rep_mean}, \code{ppp}.}
#'     \item{\code{edges}}{Data frame: \code{type}, \code{group1},
#'       \code{group2}, \code{i}, \code{j}, \code{observed}, \code{rep_mean},
#'       \code{rep_sd}, \code{z} (\code{(observed - rep_mean) / rep_sd}),
#'       \code{ppp} and \code{p_two} (\code{2 min(ppp, 1 - ppp)}).}
#'     \item{\code{draws}}{Data frame of observed and replicated global
#'       statistics per draw, for [plot_ppc()].}
#'     \item{\code{n_draws}, \code{use}, \code{edge_stat}, \code{alpha}}{
#'       Settings.}
#'   }
#'
#' @references
#' Gelman, A., Meng, X.-L. & Stern, H. (1996). Posterior predictive
#' assessment of model fitness via realized discrepancies. *Statistica
#' Sinica*, 6(4), 733-760.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
#'                      burnin = 200, nsave = 200, seed = 1)
#' pc <- posterior_predictive_check(fit, sim$data_list, ndraws = 100, seed = 1)
#' pc
#' plot_ppc(pc)
#'
#' @seealso [plot_ppc()], [mcmc_diagnostics()]
#' @export
posterior_predictive_check <- function(fit, data_list, ndraws = 200,
                                       use = c("draws", "mean"),
                                       edge_stat = c("auto", "pcor", "cor"),
                                       alpha = 0.05, chain = NULL, seed = NULL) {
  use <- match.arg(use)
  edge_stat <- match.arg(edge_stat)
  .check_edge_prob_fit(fit)
  chains <- .fit_chains(fit, chain)
  K <- chains[[1]]$K; p <- chains[[1]]$p
  if (!is.list(data_list) || length(data_list) != K)
    stop("data_list must be a list of the K = ", K, " data matrices used for fitting.",
         call. = FALSE)
  X <- lapply(data_list, function(x) {
    x <- as.matrix(x)
    if (ncol(x) != p) stop("Each data matrix must have p = ", p, " columns.", call. = FALSE)
    x <- scale(x, center = TRUE, scale = FALSE)
    matrix(x, nrow(x))
  })
  n <- vapply(X, nrow, integer(1))
  if (edge_stat == "auto") edge_stat <- if (all(n > p + 1)) "pcor" else "cor"
  if (edge_stat == "pcor" && any(n <= p + 1))
    stop("edge_stat = \"pcor\" needs more than p + 1 observations per group; ",
         "use edge_stat = \"cor\".", call. = FALSE)
  if (any(vapply(chains, function(ch) isTRUE(ch$prior_only), logical(1))))
    stop("fit is a prior_only fit.", call. = FALSE)

  # Precision matrices to replicate from: a list of [p, p, K] arrays
  if (use == "draws") {
    if (any(vapply(chains, function(ch) is.null(ch$C_save), logical(1))))
      stop("use = \"draws\" needs saved precision draws; refit with save_draws = TRUE ",
           "(use thin to keep about 200 draws), or use = \"mean\" for an approximate ",
           "check.", call. = FALSE)
    refs <- do.call(rbind, lapply(seq_along(chains), function(r)
      cbind(r, seq_len(dim(chains[[r]]$C_save)[4]))))
    pick <- unique(round(seq(1, nrow(refs), length.out = min(ndraws, nrow(refs)))))
    get_omega <- function(s) chains[[refs[pick[s], 1]]]$C_save[, , , refs[pick[s], 2], drop = FALSE][, , , 1]
    n_s <- length(pick)
  } else {
    C_mean <- Reduce(`+`, lapply(chains, function(ch) ch$C_mean * ch$n_post)) /
      sum(vapply(chains, function(ch) ch$n_post, numeric(1)))
    get_omega <- function(s) C_mean
    n_s <- as.integer(ndraws)
  }
  if (!is.null(seed)) set.seed(seed)

  pairs <- if (K > 1) t(utils::combn(K, 2)) else matrix(integer(0), 0, 2)
  ngp <- nrow(pairs)
  ut <- which(upper.tri(diag(p)), arr.ind = TRUE)
  npe <- nrow(ut)
  R_obs <- lapply(X, .ppc_cor, edge_stat = edge_stat)
  C_obs <- lapply(X, .ppc_cor, edge_stat = "cor")
  e_obs <- vapply(R_obs, function(R) R[ut], numeric(npe))
  d_obs <- vapply(seq_len(ngp), function(g) e_obs[, pairs[g, 1]] - e_obs[, pairs[g, 2]],
                  numeric(npe))
  acc <- function(r, c) list(sum = matrix(0, r, c), sumsq = matrix(0, r, c),
                             ge = matrix(0, r, c))
  ea <- acc(npe, K); da <- acc(npe, ngp)
  add <- function(a, rep, obs) {
    a$sum <- a$sum + rep; a$sumsq <- a$sumsq + rep^2; a$ge <- a$ge + (rep >= obs); a
  }
  rows <- vector("list", n_s)
  n_used <- 0L
  for (s in seq_len(n_s)) {
    Om <- get_omega(s)
    if (K == 1) Om <- array(Om, c(p, p, 1))
    U <- lapply(seq_len(K), function(k) tryCatch(chol(Om[, , k]), error = function(e) NULL))
    if (any(vapply(U, is.null, logical(1)))) next
    n_used <- n_used + 1L
    out <- list()
    e_rep <- matrix(0, npe, K)
    for (k in seq_len(K)) {
      Z <- matrix(stats::rnorm(p * n[k]), p, n[k])
      Xr <- t(backsolve(U[[k]], Z))
      Xr <- sweep(Xr, 2, colMeans(Xr))
      Rr <- .ppc_cor(Xr, edge_stat)
      e_rep[, k] <- Rr[ut]
      Cr <- if (edge_stat == "cor") Rr else .ppc_cor(Xr, "cor")
      gs_o <- .ppc_group_stats(C_obs[[k]]); gs_r <- .ppc_group_stats(Cr)
      rz_o <- .ppc_realized(X[[k]], U[[k]]); rz_r <- .ppc_realized(Xr, U[[k]])
      out[[k]] <- data.frame(draw = s, type = "group", group1 = k, group2 = NA_integer_,
                             statistic = c(names(gs_o), names(rz_o)),
                             observed = c(gs_o, rz_o), replicated = c(gs_r, rz_r),
                             stringsAsFactors = FALSE)
    }
    ea <- add(ea, e_rep, e_obs)
    if (ngp > 0) {
      d_rep <- vapply(seq_len(ngp), function(g) e_rep[, pairs[g, 1]] - e_rep[, pairs[g, 2]],
                      numeric(npe))
      d_rep <- matrix(d_rep, npe, ngp)
      da <- add(da, d_rep, matrix(d_obs, npe, ngp))
      out[[K + 1]] <- data.frame(draw = s, type = "pair", group1 = pairs[, 1],
                                 group2 = pairs[, 2], statistic = "sum_sq_diff",
                                 observed = colSums(matrix(d_obs, npe, ngp)^2),
                                 replicated = colSums(d_rep^2), stringsAsFactors = FALSE)
    }
    rows[[s]] <- do.call(rbind, out)
  }
  if (n_used < 2L) stop("Fewer than two usable posterior draws.", call. = FALSE)
  draws <- do.call(rbind, rows)
  rownames(draws) <- NULL

  key <- paste(draws$type, draws$group1, draws$group2, draws$statistic)
  global <- do.call(rbind, lapply(split(draws, factor(key, unique(key))), function(d)
    data.frame(type = d$type[1], group1 = d$group1[1], group2 = d$group2[1],
               statistic = d$statistic[1], observed = mean(d$observed),
               rep_mean = mean(d$replicated), ppp = mean(d$replicated >= d$observed),
               stringsAsFactors = FALSE)))
  rownames(global) <- NULL

  edge_df <- function(a, obs, type, g1, g2, col) {
    m <- a$sum[, col] / n_used
    sdv <- sqrt(pmax(a$sumsq[, col] - n_used * m^2, 0) / (n_used - 1))
    ppp <- a$ge[, col] / n_used
    data.frame(type = type, group1 = g1, group2 = g2, i = ut[, 1], j = ut[, 2],
               observed = obs[, col], rep_mean = m, rep_sd = sdv,
               z = (obs[, col] - m) / sdv, ppp = ppp, p_two = pmin(1, 2 * pmin(ppp, 1 - ppp)),
               stringsAsFactors = FALSE)
  }
  edges <- c(lapply(seq_len(K), function(k) edge_df(ea, e_obs, "group", k, NA_integer_, k)),
             lapply(seq_len(ngp), function(g)
               edge_df(da, matrix(d_obs, npe, ngp), "pair", pairs[g, 1], pairs[g, 2], g)))
  edges <- do.call(rbind, edges)
  rownames(edges) <- NULL

  structure(list(global = global, edges = edges, draws = draws, n_draws = n_used,
                 use = use, edge_stat = edge_stat, alpha = alpha),
            class = "multiggm_ppc")
}

#' @export
print.multiggm_ppc <- function(x, digits = 3, ...) {
  cat("<multiggm_ppc>", x$n_draws, if (x$use == "draws") "posterior draws"
      else "replicates from the posterior mean (approximate)",
      "| edge statistic:", x$edge_stat, "\n\n")
  g <- x$global
  g$observed <- signif(g$observed, digits); g$rep_mean <- signif(g$rep_mean, digits)
  g$ppp <- round(g$ppp, digits)
  print(g, row.names = FALSE)
  cat(sprintf("\nEdges with two-sided ppp < %g:\n", x$alpha))
  e <- x$edges
  key <- ifelse(e$type == "group", paste("group", e$group1),
                paste0("pair ", e$group1, "-", e$group2))
  tab <- tapply(e$p_two < x$alpha, factor(key, unique(key)), sum)
  tot <- tapply(e$p_two, factor(key, unique(key)), length)
  print(data.frame(set = names(tab), flagged = as.integer(tab), edges = as.integer(tot)),
        row.names = FALSE)
  cat("(ppp values concentrate near 0.5 under a correct model, so expect fewer than",
      "alpha x edges.)\n")
  invisible(x)
}

#' Plot posterior predictive checks
#'
#' @param x A \code{"multiggm_ppc"} object from
#'   [posterior_predictive_check()].
#' @param type \code{"global"} (default): for each global statistic, the
#'   distribution over draws of \eqn{T(y^{rep}) - T(y)}; the ppp is the
#'   share above zero. \code{"edges"}: histograms of the edge-level ppp
#'   values per group and pair, which should not pile up at 0 or 1.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso [posterior_predictive_check()]
#' @export
plot_ppc <- function(x, type = c("global", "edges")) {
  type <- match.arg(type)
  if (type == "global") {
    d <- x$draws
    d$panel <- ifelse(d$type == "group", sprintf("%s, group %d", d$statistic, d$group1),
                      sprintf("%s, groups %d-%d", d$statistic, d$group1, d$group2))
    d$diff <- d$replicated - d$observed
    ggplot2::ggplot(d, ggplot2::aes(x = .data$diff)) +
      ggplot2::geom_histogram(bins = 30, fill = "grey60", color = "white") +
      ggplot2::geom_vline(xintercept = 0, color = "firebrick") +
      ggplot2::facet_wrap(~ panel, scales = "free") +
      ggplot2::labs(title = "Posterior predictive check",
                    x = "T(replicated) - T(observed)", y = "Draws") +
      ggplot2::theme_minimal()
  } else {
    e <- x$edges
    e$panel <- ifelse(e$type == "group", paste("Group", e$group1),
                      paste0("Difference ", e$group1, "-", e$group2))
    ggplot2::ggplot(e, ggplot2::aes(x = .data$ppp)) +
      ggplot2::geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "grey60", color = "white") +
      ggplot2::facet_wrap(~ panel) +
      ggplot2::labs(title = paste("Edge-level posterior predictive p-values:", x$edge_stat),
                    x = "ppp", y = "Edges") +
      ggplot2::theme_minimal()
  }
}
