# Plot method for multiggm_fit objects

Produces diagnostic and summary plots for a fitted multi-GGM model.

## Usage

``` r
# S3 method for class 'multiggm_fit'
plot(
  x,
  type = c("trace_theta", "trace_edges", "pip", "network"),
  pip_threshold = 0.5,
  ...
)
```

## Arguments

- x:

  A `multiggm_fit` object returned by
  [`multiggm_mcmc`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md).

- type:

  Character string specifying the plot type:

  `"trace_theta"`

  :   Trace plots of the graph similarity parameter \\\theta\_{km}\\
      across post-burn-in iterations. Values \> 0 indicate the model is
      borrowing strength between groups. Requires K \>= 2.

  `"trace_edges"`

  :   Trace of the number of included edges per group across iterations.
      Useful for assessing convergence and stationarity of model
      complexity.

  `"pip"`

  :   Heatmap of posterior inclusion probabilities (PIP) for each edge,
      faceted by group. Uses viridis color scale.

  `"network"`

  :   Network visualization of the estimated graph at the given PIP
      threshold, using igraph. One plot per group.

- pip_threshold:

  Numeric; PIP threshold for edge selection in network plots. Default
  0.5. Ignored for other plot types.

- ...:

  Additional arguments (currently ignored).

## Value

For `"trace_theta"`, `"trace_edges"`, and `"pip"`: a `ggplot` object
(can be further customized). For `"network"`: invisible `NULL` (plots
are drawn as a side effect using igraph).

## See also

[`plot_trace()`](https://mljaniczek.github.io/multiGGMr/reference/plot_trace.md),
[`plot_pip_heatmap()`](https://mljaniczek.github.io/multiGGMr/reference/plot_pip_heatmap.md),
[`plot_network()`](https://mljaniczek.github.io/multiGGMr/reference/plot_network.md),
[`plot_roc()`](https://mljaniczek.github.io/multiGGMr/reference/plot_roc.md),
[`plot_recovery()`](https://mljaniczek.github.io/multiGGMr/reference/plot_recovery.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
plot(fit, type = "trace_edges")

plot(fit, type = "pip")

```
