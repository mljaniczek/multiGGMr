# multiGGM (development version)

## Shared tools: ggmtools

* `precision_to_pcor()`, `confusion_at_threshold()`, `roc_auc()` and
  `plot_roc()` now come from the new package ggmtools and are re-exported.
  spikeyglass re-exports the same functions, so attaching both packages no
  longer gives conflicting versions.
* **`roc_auc()` now handles tied scores correctly.** The previous version
  accumulated one edge at a time in sorted order, so tied edges (e.g. PIPs of
  exactly 0 or 1) were ordered by matrix position. That made the AUC depend
  on variable order and could bias it by several hundredths. The curve now
  steps only at distinct scores, and the AUC equals the Mann-Whitney
  statistic. With no tied scores, the AUC is unchanged.
* `roc_auc()` returns `NA` for `auc` (instead of an error) when the truth has
  no edges or no non-edges.
* `confusion_at_threshold()` returns a list (was a named numeric vector)
  with `precision` and `F1` added, and its third argument is now called
  `threshold` (was `thr`), with default 0.5. Use `cm$TPR` rather than
  `cm["TPR"]`. `plot_recovery()` returns these lists.
* `plot_roc()` now uses base graphics and returns the ROC list invisibly
  (it previously returned a ggplot object). It gains `add = TRUE` to overlay
  curves.
* `precision_to_pcor()` returns a matrix of `NA` with a warning (instead of
  an error) for a non-positive diagonal, and keeps dimnames.
* New methods for the ggmtools generics `pip_table()` (PIPs as a long edge
  table) and `pcor_estimate()` (posterior-mean partial correlations, one
  matrix per group). Scripts written with these work on spikeyglass fits
  too.

# multiGGM 0.4.0

## New features

* `multiggm_mcmc(method = "mghs")` fits the multiple graphical horseshoe
  (Busatto & Stingo 2025), adapted from github.com/cbusatto/mGHS. It
  returns a `multiggm_fit`, so `coef()`, `fitted()`, `posterior_pcor()`,
  `posterior_ci()`, `diff_prob_pcor()`, `plot_trace("theta")` and
  multi-chain fitting work. It reproduces upstream's draws exactly for the
  same seed (`inst/mghs/VALIDATION.md`). Local fixes: rejection samplers
  that could hang, stuck checks that missed `NaN`, and a cut-model indexing
  bug (`inst/mghs/UPSTREAM.md`). New options: `thin`, `save_draws = FALSE`
  (keep only running summaries), automatic restart of stuck chains.
* `diff_draws()`, `diff_summary()`, `bayes_fdr()` and `select_diff_edges()`
  for differential edges between any two groups, pooling chains by default.
  `select_diff_edges()` controls the posterior expected FDR using
  `P(|D| > delta)` (all samplers) or `P(g_k != g_m)` (G-Wishart and SSVS).

## Changes

* `pip_edges()`, `pip_diff_edge()` and `edge_counts()` stop with an
  explanatory error for mGHS fits, which have no edge indicators.
* `print()` no longer reads `adj_save` to count draws for mGHS fits.
