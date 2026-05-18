# Posterior probability of differential edge inclusion (K=2)

For two-group models, computes the posterior probability that each edge
has different inclusion status between the two groups (i.e., present in
one group but absent in the other).

## Usage

``` r
pip_diff_edge(fit, chain = 1L)
```

## Arguments

- fit:

  A `multiggm_fit` object (single chain) or `multiggm_fit_list` object.
  Must have `K = 2`.

- chain:

  Integer; which chain to use if `fit` is a `multiggm_fit_list`. Default
  1.

## Value

A symmetric numeric matrix `[p, p]` where entry `[i,j]` is the posterior
probability that edge (i,j) has different inclusion status between
groups 1 and 2. Values range from 0 to 1.

## See also

[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md),
[`diff_prob_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/diff_prob_pcor.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
pip_diff <- pip_diff_edge(fit)
sum(pip_diff[upper.tri(pip_diff)] > 0.5)
#> [1] 4
```
