# Differential-edge inference: posterior draws of group differences,
# per-edge summaries, and direct-posterior Bayesian FDR selection.

.fit_chains <- function(fit, chain = NULL) {
  if (inherits(fit, "multiggm_fit_list")) {
    if (is.null(chain)) return(fit$chains)
    return(fit$chains[chain])
  }
  if (inherits(fit, "multiggm_fit")) return(list(fit))
  stop("fit must be a multiggm_fit or multiggm_fit_list object.")
}

.check_groups <- function(groups, K) {
  groups <- as.integer(groups)
  if (length(groups) != 2L || anyNA(groups) || groups[1] == groups[2] ||
      any(groups < 1L | groups > K))
    stop(sprintf("groups must be two distinct group indices in 1..%d.", K))
  groups
}

.precision_to_pcor_slices <- function(Cs) {
  out <- Cs
  for (s in seq_len(dim(Cs)[3])) {
    d <- sqrt(pmax(diag(Cs[, , s]), .Machine$double.eps))
    R <- -Cs[, , s] / outer(d, d)
    diag(R) <- 1
    out[, , s] <- R
  }
  out
}

#' Posterior draws of the difference between two groups
#'
#' For each saved MCMC draw, computes the element-wise difference
#' \eqn{D = X_k - X_m} between groups \code{k} and \code{m}, where \eqn{X} is
#' either the partial correlation matrix or the precision matrix. Draws from
#' all chains of a \code{multiggm_fit_list} are pooled by default.
#'
#' @param fit A \code{multiggm_fit} or \code{multiggm_fit_list} object.
#' @param groups Integer vector of length 2, the groups \code{c(k, m)} to
#'   compare. Default \code{c(1, 2)}.
#' @param scale \code{"pcor"} (default) for partial correlations or
#'   \code{"precision"} for precision-matrix entries.
#' @param chain \code{NULL} (default) pools all chains; an integer vector
#'   selects chains of a \code{multiggm_fit_list}.
#'
#' @return Numeric array \code{[p, p, nsave_total]} of draws of \eqn{D}.
#'   Diagonal entries are 0.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' D <- diff_draws(fit)
#' dim(D)
#'
#' @seealso [diff_summary()], [select_diff_edges()]
#' @export
diff_draws <- function(fit, groups = c(1, 2), scale = c("pcor", "precision"),
                       chain = NULL) {
  scale <- match.arg(scale)
  chains <- .fit_chains(fit, chain)
  groups <- .check_groups(groups, chains[[1]]$K)
  draws <- lapply(chains, function(ch) {
    if (is.null(ch$C_save))
      stop("This fit did not save posterior draws (save_draws = FALSE). ",
           "Use diff_summary(), which reads the running summaries, or refit ",
           "with save_draws = TRUE.")
    Ck <- ch$C_save[, , groups[1], , drop = TRUE]
    Cm <- ch$C_save[, , groups[2], , drop = TRUE]
    if (length(dim(Ck)) == 2L) {
      Ck <- array(Ck, c(dim(Ck), 1L)); Cm <- array(Cm, c(dim(Cm), 1L))
    }
    if (scale == "pcor") {
      Ck <- .precision_to_pcor_slices(Ck)
      Cm <- .precision_to_pcor_slices(Cm)
    }
    D <- Ck - Cm
    p <- dim(D)[1]
    D[cbind(rep(seq_len(p), dim(D)[3]), rep(seq_len(p), dim(D)[3]),
            rep(seq_len(dim(D)[3]), each = p))] <- 0
    D
  })
  if (length(draws) == 1L) return(draws[[1]])
  p <- dim(draws[[1]])[1]
  array(unlist(draws, use.names = FALSE),
        dim = c(p, p, sum(vapply(draws, function(d) dim(d)[3], numeric(1)))))
}

#' Per-edge posterior summary of a group difference
#'
#' Summarizes the posterior of \eqn{D_{ij} = X_{k,ij} - X_{m,ij}} for every
#' edge \eqn{i < j}: mean, standard deviation, credible interval, the sign
#' probability and the exceedance probabilities \eqn{P(|D_{ij}| > c)}.
#'
#' @inheritParams diff_draws
#' @param deltas Numeric vector of magnitude thresholds \code{c} for the
#'   exceedance probabilities \eqn{P(|D| > c)}, on the scale given by
#'   \code{scale}.
#' @param probs Numeric vector of length 2 with the credible-interval
#'   quantiles, or \code{NULL} to skip the interval (faster).
#'
#' @details
#' \strong{Do not use \code{excludes_zero} or \code{p_sign} for FDR control.}
#' A credible interval that excludes zero is equivalent to
#' \code{p_sign} \eqn{\ge 1 - \alpha/2}, where
#' \code{p_sign} \eqn{= \max\{P(D > 0), P(D < 0)\}}. Whenever the posterior
#' of \eqn{D} is continuous, \eqn{P(D = 0) = 0}, so \code{p_sign} is not a
#' probability that the edge differs. This holds for the SSVS sampler (normal
#' spike), for the mGHS sampler (continuous shrinkage), and for G-Wishart
#' edges present in both groups. Edges with identical true values in both
#' groups then get \code{p_sign} near 1 from sampling noise alone. In a
#' simulation with mGHS (p = 100, n = 1000 per group), Bayesian FDR on
#' \code{p_sign} gave an achieved FDR of about 0.50 at a nominal 0.05. Use
#' \eqn{P(|D| > c)} with a practically meaningful \code{c} instead (see
#' [select_diff_edges()]).
#'
#' For mGHS fits with \code{save_draws = FALSE}, the summary is read from
#' running sums kept by the sampler over every post-burn-in iteration. Then
#' \code{deltas} must be a subset of the grid given at fit time, and the
#' credible-interval columns are \code{NA}.
#'
#' @return A data frame with one row per edge \eqn{i < j} and columns
#'   \code{i}, \code{j}, \code{mean}, \code{sd}, \code{lower}, \code{median},
#'   \code{upper}, \code{excludes_zero}, \code{p_sign}, and one
#'   \code{p_exceed_<c>} column per entry of \code{deltas}. Rows are in R's
#'   \code{upper.tri()} order. Attributes record \code{groups}, \code{scale},
#'   \code{deltas}, \code{probs} and \code{n_draws}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' ds <- diff_summary(fit, deltas = c(0.05, 0.1))
#' head(ds)
#'
#' @seealso [diff_draws()], [select_diff_edges()], [bayes_fdr()]
#' @export
diff_summary <- function(fit, groups = c(1, 2), scale = c("pcor", "precision"),
                         deltas = c(0.01, 0.02, 0.05, 0.10),
                         probs = c(0.025, 0.975), chain = NULL) {
  scale <- match.arg(scale)
  chains <- .fit_chains(fit, chain)
  K <- chains[[1]]$K; p <- chains[[1]]$p
  groups <- .check_groups(groups, K)
  if (!is.null(probs) && length(probs) != 2L)
    stop("probs must be NULL or a numeric vector of length 2.")
  idx <- which(upper.tri(diag(p)), arr.ind = TRUE)
  colnames(idx) <- c("i", "j")

  has_draws <- all(vapply(chains, function(ch) !is.null(ch$C_save), logical(1)))
  if (has_draws) {
    D <- diff_draws(fit, groups = groups, scale = scale, chain = chain)
    ns <- dim(D)[3]
    Dm <- matrix(D, nrow = p * p)[idx[, "i"] + p * (idx[, "j"] - 1L), , drop = FALSE]
    m <- rowMeans(Dm)
    s <- if (ns > 1L) sqrt(pmax(rowMeans(Dm^2) - m^2, 0) * ns / (ns - 1)) else rep(NA_real_, nrow(Dm))
    p_pos <- rowMeans(Dm > 0); p_neg <- rowMeans(Dm < 0)
    p_exc <- vapply(deltas, function(d) rowMeans(abs(Dm) > d), numeric(nrow(Dm)))
    q <- if (is.null(probs)) NULL else
      t(apply(Dm, 1, stats::quantile, probs = c(probs[1], 0.5, probs[2]), names = FALSE))
  } else {
    acc <- lapply(chains, .mghs_diff_acc, groups = groups, scale = scale)
    grid <- acc[[1]]$deltas
    j <- match(round(deltas, 12), round(grid, 12))
    if (anyNA(j))
      stop("Without saved draws, deltas must come from the grid set at fit time: ",
           paste(grid, collapse = ", "), ".")
    ns <- sum(vapply(acc, `[[`, numeric(1), "n"))
    sum_all <- function(nm) Reduce(`+`, lapply(acc, `[[`, nm))
    m <- sum_all("sum") / ns
    s <- sqrt(pmax(sum_all("sumsq") - ns * m^2, 0) / (ns - 1))
    p_pos <- sum_all("n_pos") / ns; p_neg <- sum_all("n_neg") / ns
    p_exc <- sum_all("n_exceed")[, j, drop = FALSE] / ns
    q <- if (is.null(probs)) NULL else matrix(NA_real_, nrow(idx), 3L)
    if (!is.null(probs))
      message("Credible intervals need saved draws; lower/median/upper are NA.")
  }
  p_exc <- matrix(p_exc, nrow = nrow(idx))

  out <- data.frame(idx, mean = m, sd = s)
  if (!is.null(q)) {
    out$lower <- q[, 1]; out$median <- q[, 2]; out$upper <- q[, 3]
    out$excludes_zero <- out$lower > 0 | out$upper < 0
  }
  out$p_sign <- pmax(p_pos, p_neg)
  for (d in seq_along(deltas)) out[[paste0("p_exceed_", deltas[d])]] <- p_exc[, d]
  rownames(out) <- NULL
  attr(out, "groups") <- groups
  attr(out, "scale") <- scale
  attr(out, "deltas") <- deltas
  attr(out, "probs") <- probs
  attr(out, "n_draws") <- ns
  out
}

#' Bayesian false discovery rate selection
#'
#' Selects the largest set of hypotheses whose posterior expected false
#' discovery rate is at most \code{alpha}, using the direct posterior
#' probability approach (Newton et al. 2004). Hypotheses are ranked by their
#' posterior probability \code{prob}; the expected FDR of the top \eqn{k} is
#' the mean of \eqn{1 - prob} over them.
#'
#' @param prob Numeric vector of posterior probabilities, or a symmetric
#'   matrix of edge probabilities (only the upper triangle is used).
#' @param alpha Target posterior expected FDR. Default 0.05.
#'
#' @details The guarantee is only as good as the calibration of \code{prob}.
#'   For differential edges, pass a probability of a real difference such as
#'   \eqn{P(|D| > c)} or \eqn{P(g_k \ne g_m)}, not a sign probability (see
#'   [diff_summary()]). Ties at the cutoff are broken by input order.
#'
#' @return A logical vector, or a symmetric logical matrix with \code{FALSE}
#'   diagonal, marking the selected hypotheses. Attributes:
#'   \code{expected_fdr} (posterior expected FDR of the selected set),
#'   \code{n_selected}, \code{cutoff} (smallest selected probability, or
#'   \code{NA}) and \code{alpha}.
#'
#' @references Newton, M.A., Noueiry, A., Sarkar, D. & Ahlquist, P. (2004).
#'   Detecting differential gene expression with a semiparametric
#'   hierarchical mixture method. *Biostatistics*, 5(2), 155-176.
#'
#' @examples
#' bayes_fdr(c(0.99, 0.97, 0.9, 0.6, 0.2), alpha = 0.05)
#'
#' @seealso [select_diff_edges()]
#' @export
bayes_fdr <- function(prob, alpha = 0.05) {
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1)
    stop("alpha must be a single number in [0, 1].")
  is_mat <- is.matrix(prob)
  if (is_mat) {
    if (nrow(prob) != ncol(prob)) stop("prob matrix must be square.")
    v <- prob[upper.tri(prob)]
  } else {
    v <- as.numeric(prob)
  }
  if (any(v < 0 | v > 1, na.rm = TRUE)) stop("prob must lie in [0, 1].")

  ord <- order(v, decreasing = TRUE, na.last = NA)
  running <- cumsum(1 - v[ord]) / seq_along(ord)
  k <- sum(cummax(running) <= alpha)
  sel <- logical(length(v))
  if (k > 0) sel[ord[seq_len(k)]] <- TRUE

  if (is_mat) {
    m <- matrix(FALSE, nrow(prob), ncol(prob), dimnames = dimnames(prob))
    m[upper.tri(m)] <- sel
    sel <- m | t(m)
  }
  attr(sel, "expected_fdr") <- if (k > 0) running[k] else 0
  attr(sel, "n_selected") <- k
  attr(sel, "cutoff") <- if (k > 0) v[ord[k]] else NA_real_
  attr(sel, "alpha") <- alpha
  sel
}

#' Select differential edges with Bayesian FDR control
#'
#' Computes a per-edge posterior probability that two groups differ and
#' selects edges with [bayes_fdr()]. Chains are pooled by default.
#'
#' @inheritParams diff_summary
#' @param alpha Target posterior expected FDR. Default 0.05.
#' @param rule \describe{
#'   \item{\code{"magnitude"}}{(default) \eqn{P(|D_{ij}| > \delta)}, where
#'     \eqn{D} is the group difference on \code{scale}. Works for every
#'     sampler. FDR is controlled for the hypothesis \eqn{|D_{ij}| > \delta}:
#'     a true difference smaller than \code{delta} counts as a null.}
#'   \item{\code{"structural"}}{\eqn{P(g_{k,ij} \ne g_{m,ij})}, the posterior
#'     probability that the edge is in exactly one of the two graphs. Needs
#'     edge-inclusion indicators, so it is not available for mGHS fits.}
#' }
#' @param delta Magnitude threshold for \code{rule = "magnitude"}, on the
#'   scale given by \code{scale}. Default 0.05.
#'
#' @return An object of class \code{"multiggm_diff_selection"}: a list with
#'   \code{selected} (symmetric logical \code{[p, p]} matrix), \code{prob}
#'   (symmetric \code{[p, p]} matrix of the edge probabilities used),
#'   \code{edges} (data frame of the selected edges with \code{i}, \code{j},
#'   \code{prob} and, for the magnitude rule, \code{mean} difference),
#'   \code{expected_fdr}, \code{n_selected}, \code{alpha}, \code{rule},
#'   \code{delta}, \code{scale} and \code{groups}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' sel <- select_diff_edges(fit, alpha = 0.1, delta = 0.1)
#' sel
#'
#' @seealso [diff_summary()], [bayes_fdr()], [pip_diff_edge()]
#' @export
select_diff_edges <- function(fit, groups = c(1, 2), alpha = 0.05,
                              rule = c("magnitude", "structural"),
                              delta = 0.05, scale = c("pcor", "precision"),
                              chain = NULL) {
  rule <- match.arg(rule)
  scale <- match.arg(scale)
  chains <- .fit_chains(fit, chain)
  p <- chains[[1]]$p
  groups <- .check_groups(groups, chains[[1]]$K)
  mean_diff <- NULL

  if (rule == "magnitude") {
    ds <- diff_summary(fit, groups = groups, scale = scale, deltas = delta,
                       probs = NULL, chain = chain)
    v <- ds[[ncol(ds)]]
    mean_diff <- ds$mean
  } else {
    for (ch in chains) .require_indicators(ch, "select_diff_edges(rule = \"structural\")")
    A <- lapply(chains, function(ch) ch$adj_save[, , groups[1], , drop = FALSE] !=
                                      ch$adj_save[, , groups[2], , drop = FALSE])
    n <- sum(vapply(A, function(a) dim(a)[4], numeric(1)))
    P <- Reduce(`+`, lapply(A, function(a) apply(a, c(1, 2), sum))) / n
    v <- P[upper.tri(P)]
    delta <- NA_real_
  }

  prob <- matrix(0, p, p)
  prob[upper.tri(prob)] <- v
  prob <- prob + t(prob)
  sel <- bayes_fdr(prob, alpha = alpha)

  idx <- which(upper.tri(prob), arr.ind = TRUE)
  keep <- sel[upper.tri(sel)]
  edges <- data.frame(i = idx[keep, 1], j = idx[keep, 2], prob = v[keep])
  if (!is.null(mean_diff)) edges$mean <- mean_diff[keep]
  edges <- edges[order(-edges$prob), , drop = FALSE]
  rownames(edges) <- NULL

  structure(list(selected = unname(matrix(as.logical(sel), p, p)),
                 prob = prob, edges = edges,
                 expected_fdr = attr(sel, "expected_fdr"),
                 n_selected = attr(sel, "n_selected"),
                 alpha = alpha, rule = rule, delta = delta,
                 scale = scale, groups = groups),
            class = "multiggm_diff_selection")
}

#' @export
print.multiggm_diff_selection <- function(x, ...) {
  cat("<multiggm_diff_selection>\n")
  cat(sprintf("  Groups %d vs %d, rule: %s", x$groups[1], x$groups[2], x$rule))
  if (x$rule == "magnitude") cat(sprintf(" (P(|D| > %g), %s scale)", x$delta, x$scale))
  cat("\n")
  cat(sprintf("  Selected %d edges at target FDR %g (posterior expected FDR %.3f)\n",
              x$n_selected, x$alpha, x$expected_fdr))
  if (x$n_selected > 0) print(utils::head(x$edges, 10))
  invisible(x)
}

.require_indicators <- function(fit, fn) {
  if (identical(fit$method, "mghs") || is.null(fit$adj_save))
    stop(fn, " needs edge-inclusion indicators, which the mGHS sampler does ",
         "not have (its shrinkage is continuous). For differential edges use ",
         "select_diff_edges(rule = \"magnitude\") or diff_summary(); for ",
         "edge strength use posterior_pcor() or posterior_ci().", call. = FALSE)
  invisible(TRUE)
}
