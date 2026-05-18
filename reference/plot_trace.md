# Trace plots for MCMC diagnostics

Creates trace plots to assess MCMC convergence. Supports plotting either
the graph similarity parameter (theta) or the number of edges per group
across saved iterations.

## Usage

``` r
plot_trace(fit, what = c("theta", "edges"))
```

## Arguments

- fit:

  A `multiggm_fit` object returned by
  [`multiggm_mcmc`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md).

- what:

  Character; what to plot:

  `"theta"`

  :   Graph similarity parameters \\\theta\_{km}\\ for each pair of
      groups. Requires K \>= 2.

  `"edges"`

  :   Number of edges per group at each saved iteration. Color-coded by
      group.

## Value

A `ggplot` object.

## See also

[`plot.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/plot.multiggm_fit.md),
[`edge_counts()`](https://mljaniczek.github.io/multiGGMr/reference/edge_counts.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
plot_trace(fit, "theta")

plot_trace(fit, "edges")

```
