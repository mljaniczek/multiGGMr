# multiGGM 0.5.0

## Posterior edge probabilities without saving draws

* Every sampler except `"ssvs_platform"` now keeps running sums over all
  post-burn-in iterations (`fit$edge_acc`) for a common family of per-edge
  posterior probabilities. Each is counted per iteration and averaged, never
  assembled from the groups' marginal summaries. They do not draw random
  numbers: chains, and every existing output, are unchanged for a given seed
  (mGHS outputs are bit-identical to 0.4.0; tested).
  - Per group: posterior mean and SD (precision and partial correlation),
    `P(|x| > c)`, and the structural probability: the PIP for G-Wishart and
    SSVS; for mGHS `P(w > 0.5)` and `E[w]` with `w = lambda^2/(1 + lambda^2)`,
    plus the cut-model `P(Z = 1)`.
  - Per pair of groups: `P(G_s != G_t)` (mGHS: `P(S_s != S_t)` and
    `P(Z_s != Z_t)`), `P(|D| > c)`, the union `P(structural difference or
    |D| > c)`, the mean and SD of `D`, and the sign probability (labelled as
    not a probability of difference).
* New `edge_probs()` (tidy data frame), `edge_prob_matrix()` (one p x p
  matrix) and `select_edges()` (per-group selection, rules `"structural"`,
  `"magnitude"`, `"cut"`, `"median_shrinkage"`).
* `select_diff_edges()` gains rules `"union"` and `"cut"`, a `decision`
  argument (`"bfdr"` or `"median"`) and `c` (alias of `delta`, whose default
  is unchanged). The structural rule now also works for mGHS (it used to
  stop with an error) and for fits without saved draws.
* `save_draws = FALSE` now works for the G-Wishart and SSVS samplers too:
  only the running sums, posterior means and the `Theta` trace are kept,
  which makes fits at p = 100-250 practical. `pip_edges()`,
  `pip_diff_edge()`, `coef()`, `fitted()`, `diff_summary()` and `summary()`
  use the running sums when draws are missing.
* `diff_deltas` (the thresholds `c`) applies to every sampler. Thresholds
  outside the grid can be computed after fitting when draws were saved.
* `bayes_fdr()` keeps tied probabilities together: the cutoff falls only
  between distinct values.
* The G-Wishart engine no longer allocates three unused draw-sized arrays.
* Documented (not changed): mGHS's structural and cut-model probabilities are
  weak for null edges, because the paper's weight `w` excludes the global
  scale `tau` and the upstream cut model draws `Z ~ Bernoulli(Phi(w - t))`
  rather than the paper's closed form. Use the magnitude rule for mGHS.

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

## Bug fixes

* SSVS and multi-platform samplers: the spike-and-slab update of the graph
  similarity `theta` (inherited from the Matlab reference) drew its
  proposal with `beta_prop` as a scale but used it as a rate in the
  acceptance ratio, and had the sign of the `(alpha - alpha_prop) * log(theta)`
  term reversed for spike-to-slab moves. The platform similarity `Phi` had
  the same sign error and mixed scale and rate in its normalizing terms.
  The ratios are now computed from explicit gamma densities (as the
  G-Wishart sampler already did). Default SSVS fits are essentially
  unchanged (`beta_prop = 1`, `w = 0.9`); fits with smaller `w`,
  `beta_prop != 1`, or the multi-platform sampler will differ.
  `beta_prop` is now a rate for every sampler.
* A partial `hyper` list no longer wipes out the other defaults for the
  G-Wishart and SSVS samplers. Unknown hyperparameter names are now an error.
* `simulate_multiggm(K = 1)` returned two groups.
* `simulate_multiggm_platform()` now shares one base graph across platforms
  (it drew an independent graph per platform); the data-generating process
  is documented, and the base and platform graphs are returned.

## Changes

* `pip_edges()`, `pip_diff_edge()`, `posterior_pcor()` and
  `posterior_precision()` now pool all chains by default (`chain = NULL`),
  like `diff_summary()` and `select_diff_edges()`; `chain = i` selects one.
* New `summary()`, `coef()`, `fitted()` and `plot()` methods for
  multi-chain fits (`multiggm_fit_list`).
* `method = "ssvs_platform"` now stops when `nchains > 1` and warns that
  `parallel` is ignored, instead of silently running one chain.
* Removed ineffective `EIGEN_USE_BLAS` / `EIGEN_USE_LAPACKE` defines from
  the mGHS sources (no change to results).
* Repository and CRAN hygiene: compiled objects untracked, old demo moved to
  `dev/`, `Depends: R (>= 4.0.0)`, `BugReports`, ASCII-only sources.
* New tests for each fix, plus reproducibility tests (same seed, parallel
  versus serial, R versus C++ engine).

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
