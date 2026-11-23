# Mask small partial correlations to zero

Posterior mean partial correlations of each group, with the edges that
[`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md)
does not select set to zero: by default the magnitude rule
\\P(\|\rho\_{ij}\| \> c)\\ with Bayesian FDR control. Use it to report
or draw each group's network.

## Usage

``` r
mask_pcor(
  fit,
  c = NULL,
  rule = "magnitude",
  alpha = 0.05,
  decision = NULL,
  level = 0.5,
  groups = NULL,
  chain = NULL
)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list`.

- c:

  Threshold on the partial-correlation scale (one value, or one per
  group), for `rule = "magnitude"`.

- rule, alpha, decision, level:

  Passed to
  [`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md).
  Any rule of
  [`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md)
  can be used, e.g. `"structural"` (G-Wishart, SSVS) or `"shrinkage"`
  (mGHS).

- groups:

  Groups to mask (default all).

- chain:

  `NULL` (default) pools chains.

## Value

An array `[p, p, K]` of masked posterior mean partial correlations (unit
diagonal), with attributes `c` (per group), `rule` and `n_selected`.

## Details

**Not for differential analysis.** Comparing masked networks turns
"selected in one group, not in the other" into a difference, even when
the two partial correlations are almost equal (one just above \\c\\, one
just below). Call differences with
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)
on the unmasked posterior, with \\c\\ from
[`calibrate_c()`](https://mljaniczek.github.io/multiGGMr/reference/calibrate_c.md).

## See also

[`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md),
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
                     burnin = 100, nsave = 100, save_draws = FALSE)
m <- mask_pcor(fit, c = 0.1)
round(m[, , 1], 2)
#>       [,1] [,2]  [,3] [,4] [,5]  [,6]
#> [1,]  1.00    0 -0.33    0 0.00  0.00
#> [2,]  0.00    1  0.00    0 0.00  0.00
#> [3,] -0.33    0  1.00    0 0.00 -0.38
#> [4,]  0.00    0  0.00    1 0.00  0.00
#> [5,]  0.00    0  0.00    0 1.00  0.39
#> [6,]  0.00    0 -0.38    0 0.39  1.00
```
