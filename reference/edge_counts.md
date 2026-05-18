# Per-iteration edge counts for convergence checking

Computes the number of edges in each group at each saved iteration.
Useful for trace plots of model complexity to assess MCMC convergence.

## Usage

``` r
edge_counts(fit, chain = 1L)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list` object.

- chain:

  Integer; which chain to use if `fit` is a `multiggm_fit_list`. Default
  1.

## Value

An integer matrix with `nsave` rows and `K` columns. Entry `[s, k]` is
the number of edges (upper triangle only) in group `k` at saved
iteration `s`. Column names are `"Group_1"`, `"Group_2"`, etc.

## See also

[`plot_trace()`](https://mljaniczek.github.io/multiGGMr/reference/plot_trace.md),
[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
ec <- edge_counts(fit)
plot(ec[,1], type = "l", ylab = "Edges", main = "Group 1 edge count trace")

```
