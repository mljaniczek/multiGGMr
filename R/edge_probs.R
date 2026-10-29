# Per-edge posterior probabilities from running sums kept by every sampler
# (src/edge_prob_acc.h), so edges can be selected without saving draws.

.default_deltas <- c(0.01, 0.02, 0.03, 0.05, 0.10)

# R version of EdgeProbAcc::add(), run over saved draws. Used by the pure-R
# engine, for c values outside the fit-time grid, and in tests.
.edge_acc_from_draws <- function(C_save, adj_save = NULL, deltas,
                                 weight_save = NULL, cut_save = NULL) {
  d <- dim(C_save); p <- d[1]; K <- d[3]; ns <- d[4]
  shrink <- !is.null(weight_save)
  ut <- which(upper.tri(diag(p)), arr.ind = TRUE)
  ii <- ut[, 1]; jj <- ut[, 2]; npe <- nrow(ut); nd <- length(deltas)
  pairs <- if (K > 1) t(utils::combn(K, 2)) else matrix(integer(0), 0, 2)
  storage.mode(pairs) <- "integer"
  ngp <- nrow(pairs)
  m2 <- function(r, c) matrix(0, r, c)
  a3 <- function(a, b, c) array(0, c(a, b, c))
  gsc <- function() list(sum = m2(npe, K), sumsq = m2(npe, K), n_exceed = a3(npe, nd, K))
  psc <- function() list(sum = m2(npe, ngp), sumsq = m2(npe, ngp), n_pos = m2(npe, ngp),
                         n_neg = m2(npe, ngp), n_exceed = a3(npe, nd, ngp),
                         n_union = a3(npe, nd, ngp))
  grp <- list(precision = gsc(), pcor = gsc(), n_struct = m2(npe, K), diag_sum = m2(p, K))
  pr <- list(precision = psc(), pcor = psc(), n_struct = m2(npe, ngp))
  if (shrink) {
    grp$w_sum <- m2(npe, K); grp$n_cut <- m2(npe, K); pr$n_cut <- m2(npe, ngp)
  }
  sc_names <- c("precision", "pcor")
  for (s in seq_len(ns)) {
    X <- list(precision = m2(npe, K), pcor = m2(npe, K))
    S <- Z <- matrix(FALSE, npe, K)
    for (k in seq_len(K)) {
      om <- C_save[, , k, s]
      dv <- diag(om); dg <- sqrt(dv)
      grp$diag_sum[, k] <- grp$diag_sum[, k] + dv
      w <- om[ut]
      X$precision[, k] <- w
      X$pcor[, k] <- -w / (dg[ii] * dg[jj])
      if (shrink) {
        ww <- weight_save[, , k, s][ut]
        grp$w_sum[, k] <- grp$w_sum[, k] + ww
        S[, k] <- ww > 0.5
        Z[, k] <- cut_save[, , k, s][ut] != 0
        grp$n_cut[, k] <- grp$n_cut[, k] + Z[, k]
      } else {
        S[, k] <- adj_save[, , k, s][ut] != 0
      }
      grp$n_struct[, k] <- grp$n_struct[, k] + S[, k]
      for (sc in sc_names) {
        x <- X[[sc]][, k]
        grp[[sc]]$sum[, k] <- grp[[sc]]$sum[, k] + x
        grp[[sc]]$sumsq[, k] <- grp[[sc]]$sumsq[, k] + x * x
        for (t in seq_len(nd))
          grp[[sc]]$n_exceed[, t, k] <- grp[[sc]]$n_exceed[, t, k] + (abs(x) > deltas[t])
      }
    }
    for (g in seq_len(ngp)) {
      h1 <- pairs[g, 1]; h2 <- pairs[g, 2]
      sd_ <- S[, h1] != S[, h2]
      pr$n_struct[, g] <- pr$n_struct[, g] + sd_
      ud <- sd_
      if (shrink) {
        cd <- Z[, h1] != Z[, h2]
        pr$n_cut[, g] <- pr$n_cut[, g] + cd
        ud <- cd
      }
      for (sc in sc_names) {
        D <- X[[sc]][, h1] - X[[sc]][, h2]
        q <- pr[[sc]]
        q$sum[, g] <- q$sum[, g] + D
        q$sumsq[, g] <- q$sumsq[, g] + D * D
        q$n_pos[, g] <- q$n_pos[, g] + (D > 0)
        q$n_neg[, g] <- q$n_neg[, g] + (D < 0)
        for (t in seq_len(nd)) {
          ex <- abs(D) > deltas[t]
          q$n_exceed[, t, g] <- q$n_exceed[, t, g] + ex
          q$n_union[, t, g] <- q$n_union[, t, g] + (ex | ud)
        }
        pr[[sc]] <- q
      }
    }
  }
  list(n = ns, p = p, K = K, deltas = as.numeric(deltas), group_pairs = pairs,
       shrinkage = shrink, group = grp, pair = pr)
}

# Posterior mean precision and partial-correlation arrays [p, p, K].
.acc_means <- function(acc) {
  p <- acc$p; K <- acc$K; n <- acc$n
  ut <- upper.tri(diag(p))
  C <- P <- array(0, c(p, p, K))
  for (k in seq_len(K)) {
    m <- matrix(0, p, p); m[ut] <- acc$group$precision$sum[, k] / n
    m <- m + t(m); diag(m) <- acc$group$diag_sum[, k] / n
    C[, , k] <- m
    m <- matrix(0, p, p); m[ut] <- acc$group$pcor$sum[, k] / n
    m <- m + t(m); diag(m) <- 1
    P[, , k] <- m
  }
  list(C_mean = C, pcor_mean = P)
}

# Sum the counts and sums of several chains' accumulators.
.pool_acc <- function(accs) {
  add <- function(a, b) if (is.list(a)) Map(add, a, b) else a + b
  out <- accs[[1]]
  for (a in accs[-1]) {
    if (!isTRUE(all.equal(a$deltas, out$deltas)))
      stop("Chains were fitted with different diff_deltas and cannot be pooled.")
    out$n <- out$n + a$n
    out$group <- add(out$group, a$group)
    out$pair <- add(out$pair, a$pair)
  }
  out
}

.acc_select_c <- function(acc, j) {
  sub <- function(x) {
    for (k in seq_along(x)) {
      if (names(x)[k] %in% c("n_exceed", "n_union")) x[[k]] <- x[[k]][, j, , drop = FALSE]
      else if (is.list(x[[k]])) x[[k]] <- sub(x[[k]])
    }
    x
  }
  acc$group <- sub(acc$group)
  acc$pair <- sub(acc$pair)
  acc$deltas <- acc$deltas[j]
  acc
}

.chain_edge_acc <- function(ch) {
  if (!is.null(ch$edge_acc)) return(ch$edge_acc)
  if (!identical(ch$method, "mghs") && !is.null(ch$C_save) && !is.null(ch$adj_save))
    return(.edge_acc_from_draws(ch$C_save, ch$adj_save, .default_deltas))
  stop("This fit has no edge-probability summaries; refit with multiGGM >= 0.5.0.",
       call. = FALSE)
}

# Pooled accumulator restricted to the thresholds c (NULL = the fit-time grid).
# Thresholds outside the grid are computed from saved draws when available.
.pooled_edge_acc <- function(chains, c = NULL) {
  accs <- lapply(chains, .chain_edge_acc)
  grid <- accs[[1]]$deltas
  if (is.null(c)) return(.pool_acc(accs))
  c <- as.numeric(c)
  j <- match(round(c, 12), round(grid, 12))
  if (!anyNA(j)) return(.acc_select_c(.pool_acc(accs), j))
  accs <- lapply(chains, function(ch) {
    shrink <- identical(ch$method, "mghs")
    if (is.null(ch$C_save) || (shrink && is.null(ch$weight_save)))
      stop("c = ", paste(c[is.na(j)], collapse = ", "), " is not in the grid set at fit ",
           "time (diff_deltas = ", paste(grid, collapse = ", "), ") and the draws ",
           "needed to compute it were not saved. Refit with diff_deltas including it.",
           call. = FALSE)
    if (shrink) .edge_acc_from_draws(ch$C_save, NULL, c, ch$weight_save, ch$cut_save)
    else .edge_acc_from_draws(ch$C_save, ch$adj_save, c)
  })
  .pool_acc(accs)
}

.check_edge_prob_fit <- function(fit) {
  if (inherits(fit, "multiggm_platform_fit"))
    stop("Edge probabilities are not yet available for method = \"ssvs_platform\".",
         call. = FALSE)
}

#' Posterior edge probabilities for every group and pair of groups
#'
#' Returns per-edge posterior summaries computed from running sums that every
#' sampler keeps over all post-burn-in iterations, so they are available
#' without saving draws (\code{save_draws = FALSE}). Each probability is
#' counted per MCMC iteration and then averaged, never assembled from the
#' groups' marginal summaries: the groups are dependent through the MRF prior
#' (G-Wishart, SSVS) or the correlation of edge scales (mGHS).
#'
#' @param fit A \code{multiggm_fit} or \code{multiggm_fit_list}.
#' @param scale \code{"precision"} (default) or \code{"pcor"}: the scale of
#'   \code{mean}, \code{sd}, \code{p_exceed_*}, \code{p_union_*} and
#'   \code{p_sign}.
#' @param c Thresholds for the magnitude and union probabilities. \code{NULL}
#'   (default) uses the grid given by \code{diff_deltas} at fit time. Values
#'   outside that grid are computed from saved draws when the fit has them
#'   (for mGHS this also needs the shrinkage-weight draws), and are an error
#'   otherwise.
#' @param chain \code{NULL} (default) pools all chains; an integer selects
#'   chains.
#'
#' @details
#' Rows with \code{type = "group"} describe group \code{group1}; rows with
#' \code{type = "pair"} describe the difference \eqn{D = x_{group1} -
#' x_{group2}} for every pair of groups. Columns:
#' \describe{
#'   \item{\code{mean}, \code{sd}}{Posterior mean and SD of the precision
#'     entry or partial correlation (group rows) or of \eqn{D} (pair rows).}
#'   \item{\code{p_struct}}{Structural probability. Group rows: the edge's
#'     posterior inclusion probability \eqn{P(G_{k,ij} = 1)} for G-Wishart and
#'     SSVS; for mGHS \eqn{P(w_{k,ij} > 0.5)}, where \eqn{w = \lambda^2 / (1 +
#'     \lambda^2)} is the shrinkage weight (the paper's \eqn{\kappa}). Pair
#'     rows: \eqn{P(G_s \ne G_t)}, or for mGHS \eqn{P(S_s \ne S_t)} with
#'     \eqn{S = 1(w > 0.5)}.}
#'   \item{\code{mean_weight}}{mGHS group rows: \eqn{E[w]}. The authors'
#'     median-probability rule selects \eqn{E[w] > 0.5} (Busatto & Stingo,
#'     Section 5).}
#'   \item{\code{p_cut}}{mGHS: probability of the cut-model indicator,
#'     \eqn{P(Z = 1)} (group rows) or \eqn{P(Z_s \ne Z_t)} (pair rows).}
#'   \item{\code{p_exceed_<c>}}{\eqn{P(|x| > c)} (group rows) or
#'     \eqn{P(|D| > c)} (pair rows).}
#'   \item{\code{p_union_<c>}}{Pair rows: \eqn{P(\textrm{structural difference
#'     or } |D| > c)}, counted per iteration. The structural part is the edge
#'     indicator for G-Wishart and SSVS and the cut-model indicator for mGHS.}
#'   \item{\code{p_sign}}{Pair rows: \eqn{\max\{P(D > 0), P(D < 0)\}}. This is
#'     \strong{not} a probability that the edge differs (see
#'     [diff_summary()]); it is reported for description only.}
#' }
#' Columns that do not apply are \code{NA}.
#'
#' \strong{mGHS structural and cut-model probabilities are weak.} Following
#' the paper, \eqn{w = \lambda^2 / (1 + \lambda^2)} uses the local scale only.
#' The horseshoe shrinks null edges mainly through the global scale
#' \eqn{\tau}, so for a null edge \eqn{\lambda} stays close to its
#' half-Cauchy prior, under which \eqn{P(w > 0.5) = 0.5}. Hence
#' \code{p_struct} and \code{mean_weight} are near 0.5 for null edges, and
#' \eqn{P(S_s \ne S_t)} is near 0.5 for edges that are null in both groups.
#' The cut-model indicator is drawn as \eqn{Z \sim
#' \textrm{Bernoulli}(\Phi(w - t))} (the upstream code; the paper uses the
#' closed-form \eqn{P(w > t \mid \cdot)}), so \eqn{P(Z = 1)} lies in about
#' (0.16, 0.84) and \eqn{P(Z_s \ne Z_t)} is near 0.5 for most edges. For
#' differential edges with mGHS, use the magnitude rule. The authors select
#' edges by \eqn{E[w] > 0.5} (\code{rule = "median_shrinkage"}) or by
#' \eqn{E[w] > \hat t_k}, with \eqn{\hat t_k} the posterior mean threshold
#' (attribute \code{t_hat}).
#'
#' For mGHS the result also carries the posterior mean cut-model threshold
#' per group as attribute \code{t_hat}, for the authors' rule
#' \eqn{E[w] > \hat t_k}.
#'
#' @return A data frame with columns \code{type}, \code{group1},
#'   \code{group2}, \code{i}, \code{j} and the quantities above. Edges are in
#'   R's \code{upper.tri()} order. Attributes: \code{scale}, \code{c},
#'   \code{n_draws}, \code{method} and, for mGHS, \code{t_hat}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
#'                      burnin = 100, nsave = 100, save_draws = FALSE)
#' ep <- edge_probs(fit)
#' head(ep[ep$type == "pair", ])
#'
#' @seealso [edge_prob_matrix()], [select_edges()], [select_diff_edges()]
#' @export
edge_probs <- function(fit, scale = c("precision", "pcor"), c = NULL, chain = NULL) {
  scale <- match.arg(scale)
  .check_edge_prob_fit(fit)
  chains <- .fit_chains(fit, chain)
  acc <- .pooled_edge_acc(chains, c)
  n <- acc$n; p <- acc$p; K <- acc$K; shrink <- isTRUE(acc$shrinkage)
  ut <- which(upper.tri(diag(p)), arr.ind = TRUE); npe <- nrow(ut)
  cs <- acc$deltas
  sdv <- function(sum, sumsq) {
    m <- sum / n
    if (n > 1) sqrt(pmax(sumsq - n * m^2, 0) / (n - 1)) else rep(NA_real_, length(m))
  }
  rows <- list()
  g <- acc$group; gs <- g[[scale]]
  for (k in seq_len(K)) {
    df <- data.frame(type = "group", group1 = k, group2 = NA_integer_,
                     i = ut[, 1], j = ut[, 2],
                     mean = gs$sum[, k] / n, sd = sdv(gs$sum[, k], gs$sumsq[, k]),
                     p_struct = g$n_struct[, k] / n,
                     mean_weight = if (shrink) g$w_sum[, k] / n else NA_real_,
                     p_cut = if (shrink) g$n_cut[, k] / n else NA_real_,
                     p_sign = NA_real_, stringsAsFactors = FALSE)
    for (t in seq_along(cs)) df[[paste0("p_exceed_", cs[t])]] <- gs$n_exceed[, t, k] / n
    for (t in seq_along(cs)) df[[paste0("p_union_", cs[t])]] <- NA_real_
    rows[[length(rows) + 1]] <- df
  }
  pr <- acc$pair; ps <- pr[[scale]]
  for (q in seq_len(nrow(acc$group_pairs))) {
    df <- data.frame(type = "pair", group1 = acc$group_pairs[q, 1],
                     group2 = acc$group_pairs[q, 2], i = ut[, 1], j = ut[, 2],
                     mean = ps$sum[, q] / n, sd = sdv(ps$sum[, q], ps$sumsq[, q]),
                     p_struct = pr$n_struct[, q] / n, mean_weight = NA_real_,
                     p_cut = if (shrink) pr$n_cut[, q] / n else NA_real_,
                     p_sign = pmax(ps$n_pos[, q], ps$n_neg[, q]) / n,
                     stringsAsFactors = FALSE)
    for (t in seq_along(cs)) df[[paste0("p_exceed_", cs[t])]] <- ps$n_exceed[, t, q] / n
    for (t in seq_along(cs)) df[[paste0("p_union_", cs[t])]] <- ps$n_union[, t, q] / n
    rows[[length(rows) + 1]] <- df
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$group1 <- as.integer(out$group1); out$group2 <- as.integer(out$group2)
  attr(out, "scale") <- scale
  attr(out, "c") <- cs
  attr(out, "n_draws") <- n
  attr(out, "method") <- chains[[1]]$method
  if (shrink) {
    ta <- do.call(cbind, lapply(chains, `[[`, "ta_save"))
    attr(out, "t_hat") <- rowMeans(ta)
  }
  out
}

#' One posterior edge probability as a p x p matrix
#'
#' Extracts a single quantity from [edge_probs()] for one group or one pair
#' of groups, as a symmetric matrix.
#'
#' @inheritParams edge_probs
#' @param quantity One of \code{"p_struct"}, \code{"p_cut"},
#'   \code{"mean_weight"}, \code{"mean"}, \code{"sd"}, \code{"p_sign"},
#'   \code{"p_exceed"}, \code{"p_union"} (see [edge_probs()]).
#' @param group A single group, for per-group quantities.
#' @param pair Two groups \code{c(s, t)}, for differential quantities. For
#'   \code{quantity = "mean"} the sign follows the order (\eqn{x_s - x_t}).
#' @param c A single threshold, required for \code{"p_exceed"} and
#'   \code{"p_union"}.
#'
#' @return A symmetric \code{p x p} matrix with zero diagonal.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 100, nsave = 100)
#' edge_prob_matrix(fit, "p_union", pair = c(1, 2), c = 0.05)
#'
#' @seealso [edge_probs()], [select_edges()]
#' @export
edge_prob_matrix <- function(fit, quantity = c("p_struct", "p_cut", "mean_weight", "mean",
                                               "sd", "p_sign", "p_exceed", "p_union"),
                             group = NULL, pair = NULL, c = NULL,
                             scale = c("precision", "pcor"), chain = NULL) {
  quantity <- match.arg(quantity)
  scale <- match.arg(scale)
  if (is.null(group) == is.null(pair)) stop("Give exactly one of 'group' or 'pair'.")
  needs_c <- quantity %in% c("p_exceed", "p_union")
  if (needs_c && (is.null(c) || length(c) != 1L))
    stop("quantity = \"", quantity, "\" needs a single threshold c.")
  ep <- edge_probs(fit, scale = scale, c = if (needs_c) c else NULL, chain = chain)
  col <- if (needs_c) paste0(quantity, "_", attr(ep, "c")) else quantity
  K <- max(ep$group1)
  flip <- FALSE
  if (!is.null(group)) {
    group <- as.integer(group)
    if (length(group) != 1L || group < 1L || group > K) stop("group must be one of 1..", K, ".")
    rows <- ep$type == "group" & ep$group1 == group
  } else {
    pair <- .check_groups(pair, K)
    flip <- pair[1] > pair[2]
    rows <- ep$type == "pair" & ep$group1 == min(pair) & ep$group2 == max(pair)
  }
  v <- ep[rows, col]
  if (all(is.na(v)))
    stop("\"", quantity, "\" is not available for ",
         if (is.null(group)) "pairs of groups" else "single groups",
         " with method = \"", attr(ep, "method"), "\". See ?edge_probs.", call. = FALSE)
  if (flip && quantity == "mean") v <- -v
  p <- max(ep$j)
  m <- matrix(0, p, p)
  m[upper.tri(m)] <- v
  m + t(m)
}

# Selection core shared by select_edges() and select_diff_edges().
.select_by_prob <- function(prob, decision, alpha, is_prob = TRUE) {
  p <- nrow(prob)
  if (decision == "bfdr") {
    sel <- bayes_fdr(prob, alpha = alpha)
    efdr <- attr(sel, "expected_fdr")
    sel <- unname(matrix(as.logical(sel), p, p))
  } else {
    sel <- prob > 0.5
    diag(sel) <- FALSE
    v <- prob[upper.tri(prob)][sel[upper.tri(sel)]]
    efdr <- if (!is_prob) NA_real_ else if (length(v)) mean(1 - v) else 0
  }
  list(selected = sel, expected_fdr = efdr, n_selected = sum(sel[upper.tri(sel)]))
}

.edge_table <- function(prob, sel, extra = NULL) {
  idx <- which(upper.tri(prob), arr.ind = TRUE)
  keep <- sel[upper.tri(sel)]
  edges <- data.frame(i = idx[keep, 1], j = idx[keep, 2], prob = prob[upper.tri(prob)][keep])
  if (!is.null(extra)) edges$mean <- extra[upper.tri(extra)][keep]
  edges <- edges[order(-edges$prob), , drop = FALSE]
  rownames(edges) <- NULL
  edges
}

#' Select edges within one group by posterior probability
#'
#' Selects the edges of one group's graph from the per-edge posterior
#' probabilities of [edge_probs()], either by the median-probability rule
#' (probability > 0.5) or by Bayesian FDR control ([bayes_fdr()]).
#'
#' @inheritParams edge_probs
#' @param group The group whose graph to select. Default 1.
#' @param rule \describe{
#'   \item{\code{"structural"}}{(default) \eqn{P(G_{k,ij} = 1)} (G-Wishart,
#'     SSVS) or \eqn{P(w_{k,ij} > 0.5)} (mGHS).}
#'   \item{\code{"magnitude"}}{\eqn{P(|x_{k,ij}| > c)} on \code{scale}.}
#'   \item{\code{"cut"}}{mGHS only: \eqn{P(Z_{k,ij} = 1)}, the cut-model
#'     indicator.}
#'   \item{\code{"median_shrinkage"}}{mGHS only: the authors' median rule
#'     \eqn{E[w_{k,ij}] > 0.5}. \eqn{E[w]} is not a probability, so only
#'     \code{decision = "median"} applies and no FDR is estimated.}
#' }
#' @param c Threshold for \code{rule = "magnitude"}.
#' @param decision \code{"bfdr"} (default; \code{"median"} for
#'   \code{rule = "median_shrinkage"}) selects the largest set with posterior
#'   expected FDR at most \code{alpha}, keeping tied probabilities together;
#'   \code{"median"} selects probability > 0.5.
#' @param alpha Target posterior expected FDR for \code{decision = "bfdr"}.
#'
#' @return An object of class \code{"multiggm_edge_selection"}: a list with
#'   \code{selected} (logical \code{p x p}), \code{prob}, \code{edges},
#'   \code{expected_fdr} (posterior expected FDR of the selected set),
#'   \code{n_selected}, \code{alpha}, \code{rule}, \code{c},
#'   \code{decision}, \code{scale} and \code{group}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
#'                      burnin = 100, nsave = 100, save_draws = FALSE)
#' select_edges(fit, group = 1, decision = "median")
#'
#' @seealso [edge_probs()], [select_diff_edges()], [bayes_fdr()]
#' @export
select_edges <- function(fit, group = 1, rule = c("structural", "magnitude", "cut",
                                                  "median_shrinkage"),
                         c = NULL, decision = NULL, alpha = 0.05,
                         scale = c("precision", "pcor"), chain = NULL) {
  if (identical(rule, "union"))
    stop("rule = \"union\" compares two groups; use select_diff_edges().", call. = FALSE)
  rule <- match.arg(rule)
  scale <- match.arg(scale)
  decision <- if (is.null(decision)) {
    if (rule == "median_shrinkage") "median" else "bfdr"
  } else match.arg(decision, c("bfdr", "median"))
  is_mghs <- identical(.fit_chains(fit, chain)[[1]]$method, "mghs")
  if (rule %in% c("cut", "median_shrinkage") && !is_mghs)
    stop("rule = \"", rule, "\" needs the mGHS sampler's shrinkage weights; this fit ",
         "has edge indicators, so use rule = \"structural\".", call. = FALSE)
  if (rule == "median_shrinkage" && decision != "median")
    stop("rule = \"median_shrinkage\" uses E[w], which is not a probability; ",
         "use decision = \"median\".", call. = FALSE)
  if (rule == "magnitude" && is.null(c)) stop("rule = \"magnitude\" needs c.")
  qty <- switch(rule, structural = "p_struct", magnitude = "p_exceed",
                cut = "p_cut", median_shrinkage = "mean_weight")
  prob <- edge_prob_matrix(fit, qty, group = group, c = if (rule == "magnitude") c,
                           scale = scale, chain = chain)
  s <- .select_by_prob(prob, decision, alpha, is_prob = rule != "median_shrinkage")
  structure(list(selected = s$selected, prob = prob, edges = .edge_table(prob, s$selected),
                 expected_fdr = s$expected_fdr, n_selected = s$n_selected,
                 alpha = alpha, rule = rule, c = if (rule == "magnitude") c else NA_real_,
                 decision = decision, scale = scale, group = as.integer(group)),
            class = "multiggm_edge_selection")
}

#' @export
print.multiggm_edge_selection <- function(x, ...) {
  cat("<multiggm_edge_selection>\n")
  cat(sprintf("  Group %d, rule: %s", x$group, x$rule))
  if (x$rule == "magnitude") cat(sprintf(" (P(|x| > %g), %s scale)", x$c, x$scale))
  cat(sprintf(", decision: %s\n", x$decision))
  if (x$decision == "bfdr")
    cat(sprintf("  Selected %d edges at target FDR %g", x$n_selected, x$alpha))
  else cat(sprintf("  Selected %d edges with value > 0.5", x$n_selected))
  if (!is.na(x$expected_fdr)) cat(sprintf(" (posterior expected FDR %.3f)", x$expected_fdr))
  cat("\n")
  if (x$n_selected > 0) print(utils::head(x$edges, 10))
  invisible(x)
}
