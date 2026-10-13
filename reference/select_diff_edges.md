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
  rule = c("magnitude", "structural"),
  delta = 0.05,
  scale = c("pcor", "precision"),
  chain = NULL
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

  :   (default) \\P(\|D\_{ij}\| \> \delta)\\, where \\D\\ is the group
      difference on `scale`. Works for every sampler. FDR is controlled
      for the hypothesis \\\|D\_{ij}\| \> \delta\\: a true difference
      smaller than `delta` counts as a null.

  `"structural"`

  :   \\P(g\_{k,ij} \ne g\_{m,ij})\\, the posterior probability that the
      edge is in exactly one of the two graphs. Needs edge-inclusion
      indicators, so it is not available for mGHS fits.

- delta:

  Magnitude threshold for `rule = "magnitude"`, on the scale given by
  `scale`. Default 0.05.

- scale:

  `"pcor"` (default) for partial correlations or `"precision"` for
  precision-matrix entries.

- chain:

  `NULL` (default) pools all chains; an integer vector selects chains of
  a `multiggm_fit_list`.

## Value

An object of class `"multiggm_diff_selection"`: a list with `selected`
(symmetric logical `[p, p]` matrix), `prob` (symmetric `[p, p]` matrix
of the edge probabilities used), `edges` (data frame of the selected
edges with `i`, `j`, `prob` and, for the magnitude rule, `mean`
difference), `expected_fdr`, `n_selected`, `alpha`, `rule`, `delta`,
`scale` and `groups`.

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
#>   Groups 1 vs 2, rule: magnitude (P(|D| > 0.1), pcor scale)
#>   Selected 7 edges at target FDR 0.1 (posterior expected FDR 0.083)
#>   i j prob        mean
#> 1 1 2 1.00 -0.54583861
#> 2 2 4 1.00  0.41929256
#> 3 2 5 1.00  0.48157263
#> 4 3 4 0.98  0.49929607
#> 5 7 8 0.93  0.32794996
#> 6 4 6 0.89 -0.24818914
#> 7 6 7 0.62  0.06962387
```
