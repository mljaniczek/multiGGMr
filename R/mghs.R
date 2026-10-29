# Multiple graphical horseshoe (mGHS; Busatto & Stingo 2025).
# The sampler in src/mGHS.cpp is adapted from github.com/cbusatto/mGHS
# (commit 1c9fd70, GPL >= 2, Claudio Busatto); see inst/mghs/UPSTREAM.md.

.mghs_hyper <- function(hyper) {
  h <- .merge_hyper(hyper, .default_hyper("mghs"), "mghs")
  mu <- h$t_mean; v <- h$t_var
  if (!is.numeric(mu) || mu <= 0 || mu >= 1) stop("hyper$t_mean must lie in (0, 1).")
  if (!is.numeric(v) || v <= 0 || v >= mu * (1 - mu))
    stop("hyper$t_var must lie in (0, t_mean * (1 - t_mean)).")
  a <- mu^2 * (1 - mu) / v - mu
  h$hyp_ta <- c(a, a * (1 - mu) / mu)
  h
}

.mghs_single <- function(S_list, n_vec, burnin, nsave, thin, hyper = NULL,
                         save_draws = TRUE, diff_deltas = .default_deltas,
                         max_attempts = 3L, disp = FALSE, save_weights = FALSE, ...) {
  K <- length(S_list); p <- nrow(S_list[[1L]])
  if (K < 2L) stop("method = \"mghs\" needs at least 2 groups.")
  if (!all(vapply(S_list, function(S) all(is.finite(S)), logical(1))))
    stop("S_list contains non-finite values.")
  h <- .mghs_hyper(hyper)
  diff_deltas <- sort(unique(as.numeric(diff_deltas)))
  if (any(diff_deltas < 0)) stop("diff_deltas must be non-negative.")

  for (attempt in seq_len(max_attempts)) {
    res <- mghs_engine_cpp(
      burnin = as.integer(burnin), nsave = as.integer(nsave), thin = as.integer(thin),
      n = as.numeric(n_vec), S = lapply(S_list, as.matrix), p = as.integer(p),
      hyp_ta = h$hyp_ta, save_draws = isTRUE(save_draws),
      diff_deltas = diff_deltas, verbose = isTRUE(disp), print_every = 500L,
      save_weights = isTRUE(save_weights))
    if (res$stuck == 0L) break
    if (isTRUE(disp))
      message(sprintf("mghs attempt %d stuck at iteration %d; restarting.",
                      attempt, res$stuck_iter))
  }
  if (res$stuck != 0L) {
    why <- c("local-scale parameter alpha above 1e6 or non-finite",
             "local-scale draw non-finite",
             "global-scale parameter above 1e6 or non-finite",
             "global-scale draw non-finite",
             "local-scale parameter beta above 1e6 or non-finite")[res$stuck_reason]
    stop(sprintf("mGHS sampler got stuck in all %d attempts (last at iteration %d: %s).",
                 max_attempts, res$stuck_iter, why), call. = FALSE)
  }

  structure(list(
    K = K, p = p, method = "mghs",
    C_save = res$C_save, Sig_save = NULL, adj_save = NULL,
    R_save = res$R_save, tau_save = res$tau_save,
    C_mean = res$C_mean, pcor_mean = res$pcor_mean,
    cut_z_mean = res$cut_z_mean, ta_save = res$ta_save,
    diff_acc = res$diff_acc, edge_acc = res$edge_acc, n_post = res$n_post,
    weight_save = res$weight_save, cut_save = res$cut_save,
    attempts = attempt, reject_tries = res$reject_tries,
    save_draws = isTRUE(save_draws),
    hyper = h, call = match.call()
  ), class = "multiggm_fit")
}

# Running sums of D = X_k - X_m for one chain, oriented as groups = c(k, m).
.mghs_diff_acc <- function(fit, groups, scale) {
  da <- fit$diff_acc
  if (is.null(da) && !is.null(fit$edge_acc)) {
    ea <- fit$edge_acc
    da <- list(precision = ea$pair$precision, pcor = ea$pair$pcor,
               deltas = ea$deltas, group_pairs = ea$group_pairs)
  }
  if (is.null(da))
    stop("This fit has neither saved draws nor running difference summaries.")
  g <- which(da$group_pairs[, 1] == min(groups) & da$group_pairs[, 2] == max(groups))
  a <- da[[scale]]
  flip <- groups[1] > groups[2]
  list(n = fit$n_post, deltas = da$deltas,
       sum = (if (flip) -1 else 1) * a$sum[, g],
       sumsq = a$sumsq[, g],
       n_pos = if (flip) a$n_neg[, g] else a$n_pos[, g],
       n_neg = if (flip) a$n_pos[, g] else a$n_neg[, g],
       n_exceed = matrix(a$n_exceed[, , g], nrow = dim(a$n_exceed)[1]))
}

.summary_mghs <- function(object) {
  K <- object$K
  R_mean <- apply(object$R_save, c(1, 2), mean)
  structure(list(
    K = K, p = object$p, nsave = dim(object$R_save)[3], n_post = object$n_post,
    attempts = object$attempts, save_draws = object$save_draws,
    R_mean = R_mean, tau_mean = rowMeans(object$tau_save), hyper = object$hyper
  ), class = "summary.multiggm_mghs")
}

#' @export
print.summary.multiggm_mghs <- function(x, ...) {
  cat("multiGGM MCMC Summary (mGHS, multiple graphical horseshoe)\n")
  cat("==========================================================\n")
  cat("Groups (K):", x$K, "  |  Nodes (p):", x$p, "  |  Saved draws:", x$nsave,
      if (x$save_draws) "" else "(summaries only)", "\n")
  cat("Post-burn-in iterations:", x$n_post, "  |  Attempts:", x$attempts, "\n\n")
  cat("Posterior mean of R (correlation of edge scales across groups):\n")
  print(round(x$R_mean, 3))
  cat("\nPosterior mean global shrinkage tau per group:\n")
  print(round(x$tau_mean, 4))
  cat("\nNo edge-inclusion indicators: use select_diff_edges() or diff_summary()\n",
      "for differential edges and posterior_ci() for edge strength.\n", sep = "")
  invisible(x)
}
