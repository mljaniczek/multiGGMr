# Heatmap of posterior inclusion probabilities

Creates a heatmap of posterior inclusion probabilities (PIP) using
ggplot2 with a viridis color scale. Only the upper triangle is shown.

## Usage

``` r
plot_pip_heatmap(pip, node_labels = NULL)
```

## Arguments

- pip:

  A numeric array of PIPs: either `[p, p]` (single group) or `[p, p, K]`
  (multiple groups, as returned by
  [`pip_edges`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md)).

- node_labels:

  Optional character vector of length p with node names. If `NULL`,
  nodes are labeled 1 through p.

## Value

A `ggplot` object. For K \> 1, faceted by group.

## See also

[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md),
[`plot.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/plot.multiggm_fit.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
pip <- pip_edges(fit)
plot_pip_heatmap(pip)

```
