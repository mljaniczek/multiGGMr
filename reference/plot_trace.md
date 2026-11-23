# Trace plots for MCMC diagnostics

Trace plots of scalar summaries of the chain, overlaid by chain. All but
`"theta"` use per-iteration traces that every sampler keeps, so they
work with `save_draws = FALSE`.

## Usage

``` r
plot_trace(
  fit,
  what = c("theta", "edges", "loglik", "diff"),
  rule = NULL,
  c = NULL,
  scale = c("pcor", "precision")
)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list` object returned by
  [`multiggm_mcmc`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md);
  chains are overlaid in colour.

- what:

  Character; what to plot:

  `"theta"`

  :   Graph similarity parameters \\\theta\_{km}\\ for each pair of
      groups (mGHS: edge-scale correlations \\R\\). Requires K \>= 2.

  `"edges"`

  :   Number of edges per group (mGHS: edges with shrinkage weight \\w
      \> 0.5\\).

  `"loglik"`

  :   Gaussian log-likelihood of each group.

  `"diff"`

  :   Number of differential edges per pair of groups, by `rule`.

- rule, c, scale:

  For `what = "diff"`: which edges count as differential
  (`"structural"`, the default for G-Wishart and SSVS; `"magnitude"`,
  the default for mGHS; `"union"`; `"cut"`), the threshold `c` (default
  the largest value of `diff_deltas`) and the scale of \\D\\.

## Value

A `ggplot` object.

## See also

[`mcmc_diagnostics()`](https://mljaniczek.github.io/multiGGMr/reference/mcmc_diagnostics.md),
[`plot.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/plot.multiggm_fit.md),
[`edge_counts()`](https://mljaniczek.github.io/multiGGMr/reference/edge_counts.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
plot_trace(fit, "theta")
#> Ignoring unknown labels:
#> • colour : "Chain"

plot_trace(fit, "edges")
#> Ignoring unknown labels:
#> • colour : "Chain"

plot_trace(fit, "loglik")
#> Ignoring unknown labels:
#> • colour : "Chain"

plot_trace(fit, "diff")
#> Ignoring unknown labels:
#> • colour : "Chain"

```
