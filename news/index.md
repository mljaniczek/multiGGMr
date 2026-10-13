# Changelog

## multiGGM 0.4.0

### New features

- `multiggm_mcmc(method = "mghs")` fits the multiple graphical horseshoe
  (Busatto & Stingo 2025), adapted from github.com/cbusatto/mGHS. It
  returns a `multiggm_fit`, so
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`posterior_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md),
  [`posterior_ci()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_ci.md),
  [`diff_prob_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/diff_prob_pcor.md),
  `plot_trace("theta")` and multi-chain fitting work. It reproduces
  upstream’s draws exactly for the same seed
  (`inst/mghs/VALIDATION.md`). Local fixes: rejection samplers that
  could hang, stuck checks that missed `NaN`, and a cut-model indexing
  bug (`inst/mghs/UPSTREAM.md`). New options: `thin`,
  `save_draws = FALSE` (keep only running summaries), automatic restart
  of stuck chains.
- [`diff_draws()`](https://mljaniczek.github.io/multiGGMr/reference/diff_draws.md),
  [`diff_summary()`](https://mljaniczek.github.io/multiGGMr/reference/diff_summary.md),
  [`bayes_fdr()`](https://mljaniczek.github.io/multiGGMr/reference/bayes_fdr.md)
  and
  [`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)
  for differential edges between any two groups, pooling chains by
  default.
  [`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)
  controls the posterior expected FDR using `P(|D| > delta)` (all
  samplers) or `P(g_k != g_m)` (G-Wishart and SSVS).

### Changes

- [`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md),
  [`pip_diff_edge()`](https://mljaniczek.github.io/multiGGMr/reference/pip_diff_edge.md)
  and
  [`edge_counts()`](https://mljaniczek.github.io/multiGGMr/reference/edge_counts.md)
  stop with an explanatory error for mGHS fits, which have no edge
  indicators.
- [`print()`](https://rdrr.io/r/base/print.html) no longer reads
  `adj_save` to count draws for mGHS fits.
