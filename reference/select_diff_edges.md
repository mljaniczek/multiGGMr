# Select differential edges with Bayesian FDR control

Computes a per-edge posterior probability that two groups differ and
selects edges with
[`bayes_fdr()`](https://mljaniczek.github.io/multiGGMr/reference/bayes_fdr.md).
Chains are pooled by default.

## Usage

``` r
select_diff_edges(
  fit,
  groups = c(1, 2),
  alpha = 0.05,
  rule = c("magnitude", "structural", "union", "cut"),
  delta = 0.05,
  scale = c("pcor", "precision"),
  chain = NULL,
  decision = c("bfdr", "median"),
  c = delta
)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list` object.

- groups:

  Integer vector of length 2, the groups `c(k, m)` to compare. Default
  `c(1, 2)`.

- alpha:

  Target posterior expected FDR. Default 0.05.

- rule:

  `"magnitude"`

  :   (default) \\P(\|D\_{ij}\| \> c)\\, where \\D\\ is the group
      difference on `scale`. Works for every sampler. FDR is controlled
      for the hypothesis \\\|D\_{ij}\| \> c\\: a true difference smaller
      than `c` counts as a null.

  `"structural"`

  :   \\P(G\_{s,ij} \ne G\_{t,ij})\\, the posterior probability that the
      edge is in exactly one of the two graphs. For mGHS, \\P(S_s \ne
      S_t)\\ with \\S = 1(w \> 0.5)\\ on the shrinkage weight. Detects
      edges present in one group only, not changes in strength.

  `"union"`

  :   \\P(\textrm{structural difference or } \|D\_{ij}\| \> c)\\,
      counted per iteration (mGHS: the cut-model indicator is the
      structural part). Detects both deletions and changes in strength.

  `"cut"`

  :   mGHS only: \\P(Z_s \ne Z_t)\\ for the cut-model indicator.

  See
  [`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md)
  for the exact definitions.

- delta, c:

  Threshold \\c\\ for the magnitude and union rules, on the scale given
  by `scale`: the smallest difference worth detecting. `c` is an alias
  of `delta` (default 0.05). `c` may also be the result of
  [`calibrate_c()`](https://mljaniczek.github.io/multiGGMr/reference/calibrate_c.md),
  whose `suggested_c` is then used (pcor scale only). A fixed `c` that
  is not well above the sampling noise of \\D\\ gives many false calls;
  see
  [`calibrate_c()`](https://mljaniczek.github.io/multiGGMr/reference/calibrate_c.md).

- scale:

  `"pcor"` (default) for partial correlations or `"precision"` for
  precision-matrix entries.

- chain:

  `NULL` (default) pools all chains; an integer vector selects chains of
  a `multiggm_fit_list`.

- decision:

  `"bfdr"` (default) selects the largest set with posterior expected FDR
  at most `alpha`
  ([`bayes_fdr()`](https://mljaniczek.github.io/multiGGMr/reference/bayes_fdr.md),
  keeping tied probabilities together); `"median"` selects probability
  \> 0.5.

## Value

An object of class `"multiggm_diff_selection"`: a list with `selected`
(symmetric logical `[p, p]` matrix), `prob` (symmetric `[p, p]` matrix
of the edge probabilities used), `edges` (data frame of the selected
edges with `i`, `j`, `prob`, `q_value` (see
[`bfdr_qvalues()`](https://mljaniczek.github.io/multiGGMr/reference/bfdr_qvalues.md))
and, for the magnitude rule, `mean` difference), `expected_fdr`
(posterior expected FDR of the selected set), `n_selected`, `alpha`,
`rule`, `delta` (the threshold used, or `NA`), `scale`, `groups` and
`decision`.

## See also

[`diff_summary()`](https://mljaniczek.github.io/multiGGMr/reference/diff_summary.md),
[`bayes_fdr()`](https://mljaniczek.github.io/multiGGMr/reference/bayes_fdr.md),
[`pip_diff_edge()`](https://mljaniczek.github.io/multiGGMr/reference/pip_diff_edge.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
sel <- select_diff_edges(fit, alpha = 0.1, delta = 0.1)
sel
#> <multiggm_diff_selection>
#>   Groups 1 vs 2, rule: magnitude (c = 0.1, pcor scale), decision: bfdr
#>   Selected 6 edges at target FDR 0.1 (posterior expected FDR 0.033)
#>   i j prob    q_value       mean
#> 1 1 2 1.00 0.00000000 -0.5458386
#> 2 2 4 1.00 0.00000000  0.4192926
#> 3 2 5 1.00 0.00000000  0.4815726
#> 4 3 4 0.98 0.00500000  0.4992961
#> 5 7 8 0.93 0.01800000  0.3279500
#> 6 4 6 0.89 0.03333333 -0.2481891
```
