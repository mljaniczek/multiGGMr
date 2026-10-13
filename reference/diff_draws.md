# Posterior draws of the difference between two groups

For each saved MCMC draw, computes the element-wise difference \\D =
X_k - X_m\\ between groups `k` and `m`, where \\X\\ is either the
partial correlation matrix or the precision matrix. Draws from all
chains of a `multiggm_fit_list` are pooled by default.

## Usage

``` r
diff_draws(fit, groups = c(1, 2), scale = c("pcor", "precision"), chain = NULL)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list` object.

- groups:

  Integer vector of length 2, the groups `c(k, m)` to compare. Default
  `c(1, 2)`.

- scale:

  `"pcor"` (default) for partial correlations or `"precision"` for
  precision-matrix entries.

- chain:

  `NULL` (default) pools all chains; an integer vector selects chains of
  a `multiggm_fit_list`.

## Value

Numeric array `[p, p, nsave_total]` of draws of \\D\\. Diagonal entries
are 0.

## See also

[`diff_summary()`](https://mljaniczek.github.io/multiGGMr/reference/diff_summary.md),
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
D <- diff_draws(fit)
dim(D)
#> [1]   8   8 100
```
