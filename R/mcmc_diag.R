# MCMC diagnostics, node-level and network-level differential summaries,
# WAIC, prior-versus-posterior edge probabilities and FDR calibration. All
# are computed from running sums and small per-iteration traces kept by the
# samplers (src/edge_prob_acc.h), so they work with save_draws = FALSE.

# Options passed to the C++ engines for the diagnostics.
.diag_opts <- function(S_list, n_vec, X_list, waic, prior_only, mcse_batches) {
  mcse_batches <- suppressWarnings(as.integer(mcse_batches))
  if (length(mcse_batches) != 1L || is.na(mcse_batches) || mcse_batches < 0L)
    stop("mcse_batches must be a single non-negative integer.", call. = FALSE)
  if (isTRUE(prior_only)) return(list(n_batch = mcse_batches))
  if (isTRUE(waic) && is.null(X_list))
    stop("waic = TRUE needs the observations: give data_list, not S_list.", call. = FALSE)
  dbl <- function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m }
  list(n_batch = mcse_batches, S = lapply(S_list, dbl), n = as.numeric(n_vec),
       X = if (isTRUE(waic)) lapply(X_list, dbl) else NULL)
}

# ---------------------------------------------------------------------------
# Rank-normalized split R-hat and effective sample size (Vehtari, Gelman,
# Simpson, Carpenter & Buerkner 2021), for a matrix of draws [iterations,
# chains]. Same algorithm as the posterior package.

.split_chains <- function(x) {
  n <- nrow(x)
  if (n < 4L) return(x)
  h <- n %/% 2L
  cbind(x[seq_len(h), , drop = FALSE], x[(n - h + 1L):n, , drop = FALSE])
}

.z_scale <- function(x) {
  r <- rank(x, ties.method = "average")
  z <- stats::qnorm((r - 3 / 8) / (length(x) - 2 * 3 / 8 + 1))
  matrix(z, nrow(x), ncol(x))
}

.is_constant <- function(x) {
  x <- x[is.finite(x)]
  length(x) == 0L || diff(range(x)) < sqrt(.Machine$double.eps) * max(1, abs(x))
}

.rhat_basic <- function(x) {
  if (anyNA(x) || .is_constant(x) || nrow(x) < 2L) return(NA_real_)
  n <- nrow(x)
  B <- n * stats::var(colMeans(x))
  W <- mean(apply(x, 2, stats::var))
  sqrt(((n - 1) / n * W + B / n) / W)
}

.autocov <- function(x) {
  n <- length(x)
  m2 <- 2L * stats::nextn(n)
  y <- c(x - mean(x), rep.int(0, m2 - n))
  f <- stats::fft(y)
  Re(stats::fft(Conj(f) * f, inverse = TRUE))[seq_len(n)] / (m2 * n)
}

.ess_basic <- function(x) {
  n <- nrow(x); m <- ncol(x)
  if (n < 4L || anyNA(x) || .is_constant(x)) return(NA_real_)
  acov <- apply(x, 2, .autocov)
  mean_var <- mean(acov[1, ]) * n / (n - 1)
  var_plus <- mean_var * (n - 1) / n
  if (m > 1L) var_plus <- var_plus + stats::var(colMeans(x))
  rho <- numeric(n)
  t <- 0L
  rho_even <- 1
  rho[1] <- rho_even
  rho_odd <- 1 - (mean_var - mean(acov[2, ])) / var_plus
  rho[2] <- rho_odd
  while (t < n - 5L && !is.nan(rho_even + rho_odd) && rho_even + rho_odd > 0) {
    t <- t + 2L
    rho_even <- 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
    rho_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
    if (rho_even + rho_odd >= 0) {
      rho[t + 1] <- rho_even
      rho[t + 2] <- rho_odd
    }
  }
  max_t <- t
  if (rho_even > 0) rho[max_t + 1] <- rho_even
  # Geyer's initial monotone sequence
  t <- 0L
  while (t <= max_t - 4L) {
    t <- t + 2L
    if (rho[t + 1] + rho[t + 2] > rho[t - 1] + rho[t]) {
      rho[t + 1] <- (rho[t - 1] + rho[t]) / 2
      rho[t + 2] <- rho[t + 1]
    }
  }
  ess <- m * n
  tau <- -1 + 2 * sum(rho[seq_len(max_t)]) + rho[max_t + 1]
  tau <- max(tau, 1 / log10(ess))
  ess / tau
}

# R-hat, bulk and tail ESS and MCSE of the mean for x [iterations, chains].
.convergence <- function(x) {
  x <- as.matrix(x)
  if (.is_constant(x))
    return(c(mean = mean(x), sd = 0, rhat = NA_real_, ess_bulk = NA_real_,
             ess_tail = NA_real_, mcse_mean = 0))
  xs <- .split_chains(x)
  fold <- abs(x - stats::median(x))
  rhat <- max(.rhat_basic(.z_scale(xs)), .rhat_basic(.z_scale(.split_chains(fold))),
              na.rm = TRUE)
  ess_q <- function(prob) .ess_basic(.split_chains((x <= stats::quantile(x, prob)) * 1))
  ess_mean <- .ess_basic(xs)
  c(mean = mean(x), sd = stats::sd(x), rhat = if (is.finite(rhat)) rhat else NA_real_,
    ess_bulk = .ess_basic(.z_scale(xs)),
    ess_tail = suppressWarnings(min(ess_q(0.05), ess_q(0.95), na.rm = TRUE)),
    mcse_mean = stats::sd(x) / sqrt(ess_mean))
}

# Column bind traces of equal length from several chains; shorter chains are
# cut to the shortest.
.bind_chains <- function(vs) {
  n <- min(lengths(vs))
  do.call(cbind, lapply(vs, function(v) v[seq_len(n)]))
}

.require_trace <- function(chains) {
  if (any(vapply(chains, function(ch) is.null(ch$trace), logical(1))))
    stop("This fit has no per-iteration traces; refit with multiGGM >= 0.6.0.",
         call. = FALSE)
}

.diff_trace <- function(tr, rule, t, g, scale) {
  switch(rule,
         structural = tr$n_diff_struct[, g],
         cut = {
           if (is.null(tr$n_diff_cut)) stop("rule = \"cut\" needs an mGHS fit.", call. = FALSE)
           tr$n_diff_cut[, g]
         },
         magnitude = tr$n_diff_exceed[[scale]][, t, g],
         union = tr$n_diff_union[[scale]][, t, g])
}

.grid_index <- function(acc, c) {
  j <- match(round(c, 12), round(acc$deltas, 12))
  if (is.na(j))
    stop("c = ", c, " is not in the grid set at fit time (diff_deltas = ",
         paste(acc$deltas, collapse = ", "), ").", call. = FALSE)
  j
}

.pair_index <- function(acc, groups) {
  groups <- .check_groups(groups, acc$K)
  g <- which(acc$group_pairs[, 1] == min(groups) & acc$group_pairs[, 2] == max(groups))
  g
}

.check_diff_rule <- function(rule, c, shrink) {
  if (rule %in% c("magnitude", "union") && (is.null(c) || length(c) != 1L))
    stop("rule = \"", rule, "\" needs a single threshold c.", call. = FALSE)
  if (rule == "cut" && !shrink)
    stop("rule = \"cut\" needs the mGHS sampler's cut-model indicators.", call. = FALSE)
}

#' MCMC convergence diagnostics
#'
#' Collects convergence diagnostics for scalar summaries of the chain, the
#' Monte Carlo error of edge probabilities, agreement between chains and the
#' samplers' acceptance rates in one table. Works with
#' \code{save_draws = FALSE}: the traces and batch counts it uses are always
#' kept.
#'
#' @param fit A \code{multiggm_fit} or \code{multiggm_fit_list}. Several
#'   chains are strongly recommended: R-hat from a single chain only compares
#'   its two halves.
#' @param c Threshold for the trace of the number of differential edges
#'   \eqn{|D| > c} and for the Monte Carlo error of \eqn{P(|D| > c)}.
#'   Default: the largest value of \code{diff_deltas}.
#' @param scale Scale of \eqn{D}: \code{"pcor"} (default) or
#'   \code{"precision"}.
#' @param rhat_max,ess_min,mcse_max Thresholds for the flags: R-hat above
#'   \code{rhat_max} (default 1.01), bulk or tail ESS below \code{ess_min}
#'   (default 400) and edge-probability Monte Carlo SE above \code{mcse_max}
#'   (default 0.05).
#'
#' @details
#' \strong{Scalar traces.} For each chain the samplers record, at every
#' post-burn-in iteration, the number of edges per group, the number of
#' differential edges per pair of groups (structural, \eqn{|D| > c} and union
#' for every \code{c} in \code{diff_deltas}) and the Gaussian log-likelihood
#' of each group. These, together with \eqn{\theta} (G-Wishart, SSVS) or the
#' edge-scale correlations \eqn{R} and global scales \eqn{\tau} (mGHS), are
#' summarized by rank-normalized split R-hat and bulk and tail effective
#' sample sizes (Vehtari et al. 2021), the same algorithm as the
#' \pkg{posterior} package. A whole-graph summary such as the number of
#' edges is a sensitive check: chains stuck in different parts of graph
#' space disagree on it.
#'
#' \strong{Edge probabilities.} R-hat is not informative for single 0/1 edge
#' indicators. Instead, \code{edges} reports the batch-means Monte Carlo SE
#' of every edge probability used for selection (largest and median over
#' edges, and how many exceed \code{mcse_max}), and, for several chains, the
#' largest difference in an edge probability between any two chains. An edge
#' whose probability is within a few Monte Carlo SEs of a selection cutoff is
#' not reliably selected or rejected; run the chain longer.
#'
#' @return An object of class \code{"multiggm_diagnostics"}, a list with
#'   \describe{
#'     \item{\code{scalars}}{Data frame with one row per scalar summary:
#'       \code{quantity}, \code{mean}, \code{sd}, \code{rhat},
#'       \code{ess_bulk}, \code{ess_tail}, \code{mcse_mean} and \code{flag}.}
#'     \item{\code{edges}}{Data frame with one row per edge probability:
#'       \code{quantity}, \code{max_mcse}, \code{median_mcse},
#'       \code{n_mcse_above}, \code{max_chain_diff}.}
#'     \item{\code{acceptance}}{Data frame of mean acceptance rates for the
#'       Metropolis-Hastings steps (G-Wishart and SSVS), with a flag outside
#'       (0.05, 0.95).}
#'     \item{\code{n_chains}, \code{n_iter}, \code{c}, \code{scale},
#'       \code{thresholds}}{Settings.}
#'   }
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. & Buerkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: an improved
#' R-hat for assessing convergence of MCMC. *Bayesian Analysis*, 16(2),
#' 667-718.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs", burnin = 200,
#'                      nsave = 200, nchains = 2, seed = 1, save_draws = FALSE)
#' mcmc_diagnostics(fit)
#'
#' @seealso [plot_trace()], [plot_chain_agreement()], [edge_probs()]
#' @export
mcmc_diagnostics <- function(fit, c = NULL, scale = c("pcor", "precision"),
                             rhat_max = 1.01, ess_min = 400, mcse_max = 0.05) {
  scale <- match.arg(scale)
  .check_edge_prob_fit(fit)
  chains <- .fit_chains(fit)
  .require_trace(chains)
  acc1 <- .chain_edge_acc(chains[[1]])
  if (is.null(c)) c <- max(acc1$deltas)
  t <- .grid_index(acc1, c)
  K <- acc1$K; pairs <- acc1$group_pairs
  method <- chains[[1]]$method
  is_mghs <- identical(method, "mghs")
  trs <- lapply(chains, `[[`, "trace")

  sc <- list()
  addq <- function(name, vs) sc[[name]] <<- .bind_chains(vs)
  for (k in seq_len(K)) {
    addq(sprintf("n_edges[%d]", k), lapply(trs, function(tr) tr$n_edges[, k]))
    if (!is.null(trs[[1]]$loglik))
      addq(sprintf("loglik[%d]", k), lapply(trs, function(tr) tr$loglik[, k]))
  }
  for (g in seq_len(nrow(pairs))) {
    lab <- sprintf("%d,%d", pairs[g, 1], pairs[g, 2])
    if (!is_mghs)
      addq(sprintf("n_diff_struct[%s]", lab), lapply(trs, function(tr) tr$n_diff_struct[, g]))
    addq(sprintf("n_diff_exceed_%g[%s]", c, lab),
         lapply(trs, function(tr) tr$n_diff_exceed[[scale]][, t, g]))
    if (!is_mghs)
      addq(sprintf("n_diff_union_%g[%s]", c, lab),
           lapply(trs, function(tr) tr$n_diff_union[[scale]][, t, g]))
    k <- pairs[g, 1]; m <- pairs[g, 2]
    if (is_mghs) addq(sprintf("R[%s]", lab), lapply(chains, function(ch) ch$R_save[k, m, ]))
    else if (!is.null(chains[[1]]$Theta_save))
      addq(sprintf("theta[%s]", lab), lapply(chains, function(ch) ch$Theta_save[k, m, ]))
  }
  if (is_mghs)
    for (k in seq_len(K)) addq(sprintf("tau[%d]", k), lapply(chains, function(ch) ch$tau_save[k, ]))

  rows <- lapply(names(sc), function(nm) {
    v <- .convergence(sc[[nm]])
    data.frame(quantity = nm, t(v), stringsAsFactors = FALSE, check.names = FALSE)
  })
  scalars <- do.call(rbind, rows)
  scalars$flag <- with(scalars, ifelse(!is.na(rhat) & rhat > rhat_max, "rhat",
                                ifelse(!is.na(ess_bulk) & pmin(ess_bulk, ess_tail, na.rm = TRUE) < ess_min,
                                       "ess", "")))
  rownames(scalars) <- NULL

  # Edge-probability Monte Carlo error and agreement between chains
  ep <- edge_probs(fit, scale = scale, c = c)
  epc <- if (length(chains) > 1L)
    lapply(seq_along(chains), function(r) edge_probs(fit, scale = scale, c = c, chain = r))
  cols <- c(p_struct = "mcse_struct", p_cut = "mcse_cut")
  cols[paste0("p_exceed_", c)] <- paste0("mcse_exceed_", c)
  cols[paste0("p_union_", c)] <- paste0("mcse_union_", c)
  edges <- list()
  for (type in c("group", "pair")) {
    for (pc in names(cols)) {
      rr <- ep$type == type
      v <- ep[rr, cols[[pc]]]
      if (all(is.na(ep[rr, pc]))) next
      cd <- if (!is.null(epc)) {
        P <- sapply(epc, function(e) e[e$type == type, pc])
        max(apply(P, 1, function(z) diff(range(z))))
      } else NA_real_
      edges[[length(edges) + 1]] <- data.frame(
        quantity = paste0(type, ": ", pc),
        max_mcse = if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE),
        median_mcse = if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE),
        n_mcse_above = sum(v > mcse_max, na.rm = TRUE),
        max_chain_diff = cd, stringsAsFactors = FALSE)
    }
  }
  edges <- do.call(rbind, edges)

  acceptance <- NULL
  if (!is_mghs && !is.null(chains[[1]]$ar_nu)) {
    ut <- function(m) m[upper.tri(m)]
    ar <- function(f) mean(vapply(chains, function(ch) mean(ut(ch[[f]])), numeric(1)))
    acceptance <- data.frame(
      step = c("gamma (between-model theta)", "theta (within-model)", "nu"),
      rate = c(if (K > 1) ar("ar_gamma") else NA_real_,
               if (K > 1) ar("ar_theta") else NA_real_, ar("ar_nu")),
      stringsAsFactors = FALSE)
    acceptance$flag <- ifelse(!is.na(acceptance$rate) &
                                (acceptance$rate < 0.05 | acceptance$rate > 0.95), "check", "")
  }

  structure(list(scalars = scalars, edges = edges, acceptance = acceptance,
                 n_chains = length(chains), n_iter = min(vapply(trs, function(tr)
                   nrow(tr$n_edges), integer(1))), c = c, scale = scale, method = method,
                 thresholds = c(rhat_max = rhat_max, ess_min = ess_min, mcse_max = mcse_max)),
            class = "multiggm_diagnostics")
}

#' @export
print.multiggm_diagnostics <- function(x, digits = 3, ...) {
  cat("<multiggm_diagnostics>", x$method, "|", x$n_chains, "chain(s) x", x$n_iter,
      "post-burn-in iterations\n")
  if (x$n_chains == 1L)
    cat("  One chain: R-hat compares its two halves only; run several chains.\n")
  cat("\nScalar summaries (flag: R-hat >", x$thresholds[["rhat_max"]], "or ESS <",
      x$thresholds[["ess_min"]], "):\n")
  s <- x$scalars
  num <- vapply(s, is.numeric, logical(1))
  s[num] <- lapply(s[num], round, digits)
  print(s, row.names = FALSE)
  cat("\nEdge probabilities (Monte Carlo SE; n_mcse_above: edges with SE >",
      x$thresholds[["mcse_max"]], "):\n")
  e <- x$edges
  num <- vapply(e, is.numeric, logical(1))
  e[num] <- lapply(e[num], round, digits)
  print(e, row.names = FALSE)
  if (!is.null(x$acceptance)) {
    cat("\nAcceptance rates:\n")
    a <- x$acceptance; a$rate <- round(a$rate, digits)
    print(a, row.names = FALSE)
  }
  n_flag <- sum(x$scalars$flag != "") + sum(x$edges$n_mcse_above > 0) +
    if (is.null(x$acceptance)) 0 else sum(x$acceptance$flag != "")
  cat(if (n_flag == 0) "\nNo flags.\n" else sprintf("\n%d flag(s); see above.\n", n_flag))
  invisible(x)
}

#' Plot edge probabilities of one chain against another
#'
#' Scatter plot of per-edge posterior probabilities from each chain against
#' chain 1. Points far from the diagonal mark edges on which the chains
#' disagree, the standard convergence check for graph samplers.
#'
#' @inheritParams edge_probs
#' @param fit A \code{multiggm_fit_list} with at least two chains.
#' @param quantity Column of [edge_probs()] without the threshold suffix:
#'   \code{"p_struct"} (default), \code{"p_exceed"}, \code{"p_union"} or
#'   \code{"p_cut"}.
#' @param type \code{"group"} (default) or \code{"pair"}.
#' @param c Threshold for \code{"p_exceed"} and \code{"p_union"}.
#'
#' @return A \code{ggplot} object, faceted by group (or pair) and chain.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs", burnin = 200,
#'                      nsave = 200, nchains = 2, seed = 1, save_draws = FALSE)
#' plot_chain_agreement(fit)
#'
#' @seealso [mcmc_diagnostics()]
#' @export
plot_chain_agreement <- function(fit, quantity = c("p_struct", "p_exceed", "p_union", "p_cut"),
                                 type = c("group", "pair"), c = NULL,
                                 scale = c("pcor", "precision")) {
  quantity <- match.arg(quantity)
  type <- match.arg(type)
  scale <- match.arg(scale)
  chains <- .fit_chains(fit)
  if (length(chains) < 2L) stop("plot_chain_agreement() needs at least two chains.")
  needs_c <- quantity %in% c("p_exceed", "p_union")
  if (needs_c && (is.null(c) || length(c) != 1L))
    stop("quantity = \"", quantity, "\" needs a single threshold c.")
  col <- if (needs_c) paste0(quantity, "_", c) else quantity
  eps <- lapply(seq_along(chains), function(r)
    edge_probs(fit, scale = scale, c = if (needs_c) c, chain = r))
  if (!col %in% names(eps[[1]])) col <- paste0(quantity, "_", attr(eps[[1]], "c"))
  ref <- eps[[1]][eps[[1]]$type == type, ]
  lab <- if (type == "group") paste("Group", ref$group1)
         else paste0("Groups ", ref$group1, "-", ref$group2)
  df <- do.call(rbind, lapply(seq_along(eps)[-1], function(r) {
    e <- eps[[r]][eps[[r]]$type == type, ]
    data.frame(chain1 = ref[[col]], other = e[[col]], panel = lab,
               chain = paste("Chain", r), stringsAsFactors = FALSE)
  }))
  if (all(is.na(df$chain1)))
    stop("\"", quantity, "\" is not available for this fit and type.", call. = FALSE)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$chain1, y = .data$other)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggplot2::facet_grid(chain ~ panel) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = paste("Chain agreement:", col), x = "Chain 1", y = "Other chain") +
    ggplot2::theme_minimal()
}

#' Node-level differential network summaries
#'
#' For each node, posterior summaries of how much its connections differ
#' between two groups, counted per MCMC iteration: the expected number of
#' differential edges at the node (its differential degree), the probability
#' that it has at least one, and its expected rewiring score
#' \eqn{\sum_j |D_{ij}|}. These identify hub changes that per-edge
#' selection can miss when many small changes meet at one node.
#'
#' @inheritParams select_diff_edges
#' @param rule Which edges count as differential: \code{"magnitude"}
#'   (\eqn{|D| > c}, default), \code{"union"} (structural difference or
#'   \eqn{|D| > c}), \code{"structural"} or \code{"cut"} (mGHS).
#' @param c Threshold for \code{"magnitude"} and \code{"union"}; must be in
#'   the fit's \code{diff_deltas} grid unless draws were saved.
#' @param scale Scale of \eqn{D} and of the rewiring score: \code{"pcor"}
#'   (default) or \code{"precision"}.
#'
#' @return A data frame with one row per node: \code{node},
#'   \code{diff_degree} (posterior mean number of differential edges),
#'   \code{p_any} (posterior probability of at least one), \code{rewiring}
#'   (posterior mean of \eqn{\sum_j |D_{ij}|} on \code{scale}) and
#'   \code{rewiring_rank}. Attributes \code{groups}, \code{rule}, \code{c}
#'   and \code{scale}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
#'                      burnin = 100, nsave = 100, save_draws = FALSE)
#' node_diff_summary(fit, c = 0.1)
#'
#' @seealso [diff_count_posterior()], [select_diff_edges()]
#' @export
node_diff_summary <- function(fit, groups = c(1, 2),
                              rule = c("magnitude", "union", "structural", "cut"),
                              c = NULL, scale = c("pcor", "precision"), chain = NULL) {
  rule <- match.arg(rule)
  scale <- match.arg(scale)
  .check_edge_prob_fit(fit)
  chains <- .fit_chains(fit, chain)
  acc <- .pooled_edge_acc(chains, if (rule %in% c("magnitude", "union")) c)
  .check_diff_rule(rule, c, isTRUE(acc$shrinkage))
  if (is.null(acc$node))
    stop("This fit has no node-level summaries; refit with multiGGM >= 0.6.0.", call. = FALSE)
  g <- .pair_index(acc, groups)
  nd <- acc$node; n <- acc$n
  deg <- switch(rule,
                structural = list(nd$deg_struct[, g], nd$n_any_struct[, g]),
                cut = list(nd$deg_cut[, g], nd$n_any_cut[, g]),
                magnitude = list(nd[[scale]]$deg_exceed[, 1, g], nd[[scale]]$n_any_exceed[, 1, g]),
                union = list(nd[[scale]]$deg_union[, 1, g], nd[[scale]]$n_any_union[, 1, g]))
  rw <- nd[[scale]]$rewire_sum[, g] / n
  out <- data.frame(node = seq_len(acc$p), diff_degree = deg[[1]] / n, p_any = deg[[2]] / n,
                    rewiring = rw, rewiring_rank = rank(-rw, ties.method = "min"))
  attr(out, "groups") <- as.integer(groups)
  attr(out, "rule") <- rule
  attr(out, "c") <- if (rule %in% c("magnitude", "union")) c else NA_real_
  attr(out, "scale") <- scale
  out
}

#' Posterior distribution of the number of differential edges
#'
#' The number of edges that differ between two groups, recorded at every
#' post-burn-in iteration. Its posterior gives a network-level answer to
#' "do these networks differ?" (\code{p_any}, the posterior probability of at
#' least one differential edge) and an estimate with uncertainty of how much
#' they differ, without choosing an FDR level.
#'
#' @inheritParams node_diff_summary
#' @param c Threshold for \code{"magnitude"} and \code{"union"}.
#'
#' @details
#' With continuous priors (SSVS spike, mGHS) every \eqn{D} is nonzero, so the
#' count is meaningful only with a threshold \code{c} above sampling noise
#' (see the edge-selection vignette). Compare the posterior with the same
#' count under the prior (fit with \code{prior_only = TRUE}). The posterior
#' odds of \code{p_any} divided by its prior odds is a network-level Bayes
#' factor for any difference, informative only when the prior probability
#' is not close to 1 (with many edges and a small \code{c} the prior almost
#' surely has some difference).
#'
#' @return An object of class \code{"multiggm_diff_count"}: a list with
#'   \code{draws} (the count at every iteration, chains concatenated),
#'   \code{mean}, \code{sd}, \code{quantiles} (2.5\%, 50\%, 97.5\%),
#'   \code{p_any}, \code{groups}, \code{rule}, \code{c} and \code{scale}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
#'                      burnin = 100, nsave = 100, save_draws = FALSE)
#' diff_count_posterior(fit, rule = "union", c = 0.1)
#'
#' @seealso [node_diff_summary()], [select_diff_edges()]
#' @export
diff_count_posterior <- function(fit, groups = c(1, 2),
                                 rule = c("magnitude", "union", "structural", "cut"),
                                 c = NULL, scale = c("pcor", "precision"), chain = NULL) {
  rule <- match.arg(rule)
  scale <- match.arg(scale)
  .check_edge_prob_fit(fit)
  chains <- .fit_chains(fit, chain)
  .require_trace(chains)
  acc1 <- .chain_edge_acc(chains[[1]])
  .check_diff_rule(rule, c, isTRUE(acc1$shrinkage))
  g <- .pair_index(acc1, groups)
  draws <- unlist(lapply(chains, function(ch) {
    tr <- ch$trace; t <- 1L
    if (rule %in% c("magnitude", "union")) {
      t <- match(round(c, 12), round(ch$edge_acc$deltas, 12))
      if (is.na(t)) {
        shrink <- identical(ch$method, "mghs")
        if (is.null(ch$C_save) || (shrink && is.null(ch$weight_save)))
          stop("c = ", c, " is not in the grid set at fit time (diff_deltas = ",
               paste(ch$edge_acc$deltas, collapse = ", "), ") and the draws needed ",
               "to compute it were not saved.", call. = FALSE)
        tr <- .diag_from_draws(ch$C_save, ch$adj_save, c, ch$weight_save, ch$cut_save)$trace
        t <- 1L
      }
    }
    .diff_trace(tr, rule, t, g, scale)
  }))
  structure(list(draws = draws, mean = mean(draws), sd = stats::sd(draws),
                 quantiles = stats::quantile(draws, c(0.025, 0.5, 0.975)),
                 p_any = mean(draws > 0), groups = as.integer(groups), rule = rule,
                 c = if (rule %in% c("magnitude", "union")) c else NA_real_, scale = scale),
            class = "multiggm_diff_count")
}

#' @export
print.multiggm_diff_count <- function(x, ...) {
  cat("<multiggm_diff_count> groups", x$groups[1], "vs", x$groups[2], "| rule:", x$rule)
  if (!is.na(x$c)) cat(sprintf(" (c = %g, %s scale)", x$c, x$scale))
  cat(sprintf("\n  Differential edges: mean %.1f, sd %.1f, 95%% interval [%g, %g]\n",
              x$mean, x$sd, x$quantiles[1], x$quantiles[3]))
  cat(sprintf("  P(at least one differential edge) = %.3f\n", x$p_any))
  invisible(x)
}

#' Bayesian FDR q-values
#'
#' For each edge, the smallest target \code{alpha} at which [bayes_fdr()]
#' would select it: the posterior expected FDR of the set of edges with
#' probability at least as high (ties kept together). Selecting the edges
#' with \code{q <= alpha} gives exactly the selection of
#' \code{bayes_fdr(prob, alpha)}.
#'
#' @param prob A vector of posterior probabilities, or a square matrix (the
#'   upper triangle is used).
#'
#' @return Values of the same shape as \code{prob} (a symmetric matrix with
#'   \code{NA} diagonal for matrix input).
#'
#' @examples
#' bfdr_qvalues(c(0.99, 0.95, 0.6, 0.2))
#'
#' @seealso [bayes_fdr()]
#' @export
bfdr_qvalues <- function(prob) {
  is_mat <- is.matrix(prob)
  v <- if (is_mat) prob[upper.tri(prob)] else as.numeric(prob)
  if (any(v < 0 | v > 1, na.rm = TRUE)) stop("prob must lie in [0, 1].")
  q <- rep(NA_real_, length(v))
  ord <- order(v, decreasing = TRUE, na.last = NA)
  if (length(ord)) {
    vs <- v[ord]
    running <- cummax(cumsum(1 - vs) / seq_along(vs))
    ends <- c(which(diff(vs) != 0), length(vs))
    grp <- findInterval(seq_along(vs), c(0, ends) + 1, left.open = FALSE)
    q[ord] <- running[ends[grp]]
  }
  if (is_mat) {
    m <- matrix(NA_real_, nrow(prob), ncol(prob), dimnames = dimnames(prob))
    m[upper.tri(m)] <- q
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
  q
}

#' WAIC from running sums
#'
#' The widely applicable information criterion (Watanabe 2010) of the
#' Gaussian likelihood, computed from per-observation running sums kept
#' during sampling, so it needs no saved draws. Fit with
#' \code{waic = TRUE} and \code{data_list}.
#'
#' @param fit A \code{multiggm_fit} or \code{multiggm_fit_list} fitted with
#'   \code{waic = TRUE}.
#' @param chain \code{NULL} (default) pools all chains.
#'
#' @details
#' For observation \eqn{i} with log-likelihood \eqn{\ell_{is}} at iteration
#' \eqn{s = 1, \dots, S}: \eqn{\textrm{lppd}_i = \log(S^{-1} \sum_s
#' e^{\ell_{is}})}, \eqn{p_i = \textrm{var}_s(\ell_{is})} and
#' \eqn{\textrm{elpd}_i = \textrm{lppd}_i - p_i}. WAIC is \eqn{-2 \sum_i
#' \textrm{elpd}_i}, lower is better. Use it to compare fits to the same
#' \code{data_list}, for example different \code{hyper} settings, samplers,
#' or a joint fit against separate fits per group (compare the per-group
#' rows). The data are column-centered as for fitting. Observations with
#' \eqn{p_i > 0.4} make WAIC unreliable (a warning is given); PSIS-LOO from
#' the \pkg{loo} package is then preferable and needs saved draws.
#'
#' @return An object of class \code{"multiggm_waic"}: a list with
#'   \code{estimates} (data frame with rows \code{elpd_waic}, \code{p_waic},
#'   \code{waic} and columns \code{estimate}, \code{se}), \code{by_group}
#'   (the same per group) and \code{pointwise} (data frame with
#'   \code{group}, \code{obs}, \code{elpd_waic}, \code{p_waic}).
#'
#' @references
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' *Journal of Machine Learning Research*, 11, 3571-3594.
#'
#' Vehtari, A., Gelman, A. & Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. *Statistics and
#' Computing*, 27, 1413-1432.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
#'                      burnin = 100, nsave = 100, save_draws = FALSE, waic = TRUE)
#' compute_waic(fit)
#'
#' @export
compute_waic <- function(fit, chain = NULL) {
  chains <- .fit_chains(fit, chain)
  ws <- lapply(chains, `[[`, "waic_acc")
  if (any(vapply(ws, is.null, logical(1))))
    stop("This fit has no WAIC sums; refit with waic = TRUE and data_list.", call. = FALSE)
  K <- length(ws[[1]])
  logadd <- function(a, b) ifelse(a == -Inf, b, ifelse(b == -Inf, a,
                                  pmax(a, b) + log1p(exp(-abs(a - b)))))
  pw <- lapply(seq_len(K), function(k) {
    S <- sum(vapply(ws, function(w) w[[k]]$n, numeric(1)))
    lse <- Reduce(logadd, lapply(ws, function(w) w[[k]]$lse))
    s1 <- Reduce(`+`, lapply(ws, function(w) w[[k]]$sum))
    s2 <- Reduce(`+`, lapply(ws, function(w) w[[k]]$sumsq))
    p_i <- pmax(s2 - s1^2 / S, 0) / (S - 1)
    lppd <- lse - log(S)
    data.frame(group = k, obs = seq_along(lse), elpd_waic = lppd - p_i, p_waic = p_i)
  })
  pointwise <- do.call(rbind, pw)
  est <- function(d) {
    n <- nrow(d)
    e <- sum(d$elpd_waic); pv <- sum(d$p_waic)
    se_e <- sqrt(n * stats::var(d$elpd_waic)); se_p <- sqrt(n * stats::var(d$p_waic))
    data.frame(estimate = c(e, pv, -2 * e), se = c(se_e, se_p, 2 * se_e),
               row.names = c("elpd_waic", "p_waic", "waic"))
  }
  n_bad <- sum(pointwise$p_waic > 0.4)
  if (n_bad > 0)
    warning(sprintf("%d of %d observations have p_waic > 0.4; WAIC may be unreliable.",
                    n_bad, nrow(pointwise)), call. = FALSE)
  structure(list(estimates = est(pointwise),
                 by_group = lapply(split(pointwise, pointwise$group), est),
                 pointwise = pointwise),
            class = "multiggm_waic")
}

#' @export
print.multiggm_waic <- function(x, digits = 1, ...) {
  cat("<multiggm_waic> computed from", nrow(x$pointwise), "observations\n")
  print(round(x$estimates, digits))
  invisible(x)
}

#' Prior versus posterior edge probabilities and Bayes factors
#'
#' Compares each edge probability from [edge_probs()] with the same
#' probability under the prior alone, estimated from a fit with
#' \code{prior_only = TRUE}. The prior probability of an edge, or of a
#' difference \eqn{|D| > c}, is set by the MRF and slab/spike or G-Wishart
#' priors and is generally not 0.5, so a posterior probability is only
#' interpretable next to it. The Bayes factor (posterior odds / prior odds)
#' measures how much the data moved each edge.
#'
#' @param fit The fitted model.
#' @param prior_fit A fit of the same model (same \code{method},
#'   \code{hyper}, groups, \code{p} and \code{diff_deltas}) with
#'   \code{prior_only = TRUE}.
#' @param quantity \code{"p_struct"} (default), \code{"p_exceed"} or
#'   \code{"p_union"}.
#' @param c Threshold for \code{"p_exceed"} and \code{"p_union"}.
#' @param scale Scale of \eqn{x} or \eqn{D} for the magnitude quantities.
#' @param pool_edges If \code{TRUE} (default), the prior probability is
#'   averaged over the edges of each group (or pair). The default priors
#'   treat all edges alike, so this removes Monte Carlo noise; use
#'   \code{FALSE} for edge-specific priors.
#'
#' @details
#' To keep Bayes factors finite, probabilities are estimated as
#' \eqn{(x + 0.5) / (n + 1)} from \eqn{x} of \eqn{n} iterations, which
#' bounds the Bayes factor by about \eqn{2n} times the prior odds.
#' Not available for mGHS, whose diagonal prior is improper; under its
#' half-Cauchy local scales the shrinkage weight has a Beta(1/2, 1/2) prior
#' and \eqn{P(w > 0.5) = 0.5}.
#'
#' @return A data frame with \code{type}, \code{group1}, \code{group2},
#'   \code{i}, \code{j}, \code{prior}, \code{posterior}, \code{bayes_factor}
#'   and \code{log10_bf}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
#'                      burnin = 100, nsave = 200, save_draws = FALSE, seed = 1)
#' prior <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
#'                        burnin = 100, nsave = 200, save_draws = FALSE,
#'                        seed = 2, prior_only = TRUE)
#' head(edge_bayes_factors(fit, prior))
#'
#' @export
edge_bayes_factors <- function(fit, prior_fit, quantity = c("p_struct", "p_exceed", "p_union"),
                               c = NULL, scale = c("pcor", "precision"), pool_edges = TRUE) {
  quantity <- match.arg(quantity)
  scale <- match.arg(scale)
  pch <- .fit_chains(prior_fit); fch <- .fit_chains(fit)
  if (!all(vapply(pch, function(ch) isTRUE(ch$prior_only), logical(1))))
    stop("prior_fit must be fitted with prior_only = TRUE.", call. = FALSE)
  if (pch[[1]]$method != fch[[1]]$method || pch[[1]]$p != fch[[1]]$p ||
      pch[[1]]$K != fch[[1]]$K)
    stop("fit and prior_fit must use the same method, p and number of groups.", call. = FALSE)
  if (!isTRUE(all.equal(pch[[1]]$hyper, fch[[1]]$hyper)))
    warning("fit and prior_fit have different hyperparameters.", call. = FALSE)
  needs_c <- quantity != "p_struct"
  if (needs_c && (is.null(c) || length(c) != 1L))
    stop("quantity = \"", quantity, "\" needs a single threshold c.", call. = FALSE)
  col <- if (needs_c) paste0(quantity, "_", c) else quantity
  post <- edge_probs(fit, scale = scale, c = if (needs_c) c)
  pri <- edge_probs(prior_fit, scale = scale, c = if (needs_c) c)
  col <- if (needs_c) paste0(quantity, "_", attr(post, "c")) else quantity
  n_post <- attr(post, "n_draws"); n_pri <- attr(pri, "n_draws")
  keep <- !is.na(post[[col]])
  out <- post[keep, c("type", "group1", "group2", "i", "j")]
  pr <- pri[[col]][keep]
  if (pool_edges) {
    key <- paste(out$type, out$group1, out$group2)
    pr <- stats::ave(pr, key)
  }
  adj <- function(x, n) (x * n + 0.5) / (n + 1)
  pr_a <- adj(pr, if (pool_edges) n_pri * stats::ave(pr, key, FUN = length) else n_pri)
  po_a <- adj(post[[col]][keep], n_post)
  out$prior <- pr
  out$posterior <- post[[col]][keep]
  out$bayes_factor <- (po_a / (1 - po_a)) / (pr_a / (1 - pr_a))
  out$log10_bf <- log10(out$bayes_factor)
  rownames(out) <- NULL
  attr(out, "quantity") <- col
  out
}

#' Calibration of Bayesian FDR selection against a known truth
#'
#' For simulation studies: selects edges by Bayesian FDR at each target
#' level \code{alpha} and compares the posterior expected FDR with the
#' achieved false discovery proportion and power under the true precision
#' matrices. Well-calibrated edge probabilities give achieved FDR close to
#' the expected FDR.
#'
#' @inheritParams select_diff_edges
#' @param truth A list of the K true precision matrices, or the output of
#'   [simulate_multiggm()].
#' @param group A single group to calibrate within-group selection
#'   ([select_edges()]); \code{NULL} (default) calibrates differential
#'   selection between \code{groups}.
#' @param rule Selection rule as in [select_diff_edges()] or
#'   [select_edges()].
#' @param c Threshold for the magnitude and union rules.
#' @param truth_c Threshold defining a true difference (or true edge) for
#'   the magnitude rule, on \code{scale}. Default \code{c}: a true
#'   difference is \eqn{|D_{true}| > c}. Structural truth is a difference in
#'   the nonzero pattern (or a nonzero entry); union truth is either.
#' @param alpha Vector of target FDR levels.
#' @param tol Entries with absolute value at most \code{tol} count as zero in
#'   the true precision matrices.
#'
#' @return A data frame with \code{alpha}, \code{n_selected},
#'   \code{expected_fdr} (posterior expected FDR of the selection),
#'   \code{achieved_fdr} (false discovery proportion), \code{power} and
#'   \code{n_true}. Attributes \code{rule}, \code{c}, \code{truth_c}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 100, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
#'                      burnin = 200, nsave = 200, save_draws = FALSE)
#' fdr_calibration(fit, sim, rule = "union", c = 0.1)
#'
#' @export
fdr_calibration <- function(fit, truth, groups = c(1, 2), group = NULL,
                            rule = c("magnitude", "structural", "union", "cut"),
                            c = NULL, truth_c = c, scale = c("pcor", "precision"),
                            alpha = seq(0.01, 0.20, by = 0.01), tol = 1e-8, chain = NULL) {
  rule <- match.arg(rule)
  scale <- match.arg(scale)
  if (is.list(truth) && !is.null(truth$Omega_list)) truth <- truth$Omega_list
  if (!is.list(truth)) stop("truth must be a list of true precision matrices.")
  to_scale <- function(O) if (scale == "pcor") {
    d <- sqrt(diag(O)); P <- -O / outer(d, d); diag(P) <- 1; P
  } else O
  if (rule %in% c("magnitude", "union") && (is.null(c) || length(c) != 1L))
    stop("rule = \"", rule, "\" needs a single threshold c.", call. = FALSE)
  qty <- switch(rule, structural = "p_struct", magnitude = "p_exceed",
                union = "p_union", cut = "p_cut")
  cc <- if (rule %in% c("magnitude", "union")) c
  if (is.null(group)) {
    groups <- .check_groups(groups, length(truth))
    prob <- edge_prob_matrix(fit, qty, pair = groups, c = cc, scale = scale, chain = chain)
    O1 <- truth[[groups[1]]]; O2 <- truth[[groups[2]]]
    s_true <- (abs(O1) > tol) != (abs(O2) > tol)
    m_true <- if (!is.null(truth_c)) abs(to_scale(O1) - to_scale(O2)) > truth_c
    true <- switch(rule, structural = , cut = s_true, magnitude = m_true,
                   union = s_true | m_true)
  } else {
    if (rule == "union") stop("rule = \"union\" compares two groups.", call. = FALSE)
    prob <- edge_prob_matrix(fit, qty, group = group, c = cc, scale = scale, chain = chain)
    O <- truth[[group]]
    true <- switch(rule, structural = , cut = abs(O) > tol,
                   magnitude = abs(to_scale(O)) > truth_c)
  }
  ut <- upper.tri(prob)
  tv <- true[ut]
  rows <- lapply(alpha, function(a) {
    sel <- bayes_fdr(prob, alpha = a)
    sv <- sel[ut]
    ns <- sum(sv)
    data.frame(alpha = a, n_selected = ns, expected_fdr = attr(sel, "expected_fdr"),
               achieved_fdr = if (ns > 0) sum(sv & !tv) / ns else 0,
               power = if (sum(tv) > 0) sum(sv & tv) / sum(tv) else NA_real_,
               n_true = sum(tv))
  })
  out <- do.call(rbind, rows)
  attr(out, "rule") <- rule
  attr(out, "c") <- if (is.null(cc)) NA_real_ else cc
  attr(out, "truth_c") <- if (is.null(truth_c)) NA_real_ else truth_c
  out
}
