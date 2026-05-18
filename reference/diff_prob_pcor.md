# Posterior probability of differential partial correlation (K=2)

For two-group models, computes the posterior probability that the
partial correlation between each pair of variables differs between
groups by more than a threshold `delta`.

## Usage

``` r
diff_prob_pcor(pcor, delta = 0)
```

## Arguments

- pcor:

  A numeric array of dimension `[p, p, 2, nsave]`, as returned by
  [`posterior_pcor`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md)
  with `K = 2`.

- delta:

  Numeric; practical difference threshold. Default 0 (any difference).
  Set to e.g. 0.1 to detect meaningful differences only.

## Value

A symmetric numeric matrix `[p, p]` where entry `[i,j]` is the posterior
probability that \\\|\rho\_{ij}^{(1)} - \rho\_{ij}^{(2)}\| \> \delta\\.

## See also

[`posterior_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md),
[`pip_diff_edge()`](https://mljaniczek.github.io/multiGGMr/reference/pip_diff_edge.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
pcor_draws <- posterior_pcor(fit)
diff_prob <- diff_prob_pcor(pcor_draws, delta = 0.1)
# Edges with high probability of differential partial correlation
sum(diff_prob[upper.tri(diff_prob)] > 0.5)
#> [1] 10
```
