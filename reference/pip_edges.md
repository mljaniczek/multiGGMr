# Posterior inclusion probabilities (PIP) per group

Computes the posterior inclusion probability for each edge in each
group. The PIP is the fraction of saved MCMC iterations in which each
edge was included in the graph.

## Usage

``` r
pip_edges(fit, chain = NULL, platform = NULL)
```

## Arguments

- fit:

  A `multiggm_fit` object (single chain), `multiggm_fit_list` object, or
  `multiggm_platform_fit` object.

- chain:

  `NULL` (default) pools the draws of all chains of a
  `multiggm_fit_list`; an integer vector selects chains.

- platform:

  Integer; which platform to extract (required for
  `multiggm_platform_fit` objects). Ignored for single-platform fits.

## Value

A numeric array of dimension `[p, p, K]`. Entry `pip[i, j, k]` is the
posterior probability that edge (i,j) is present in the graph of group
`k`. Values range from 0 to 1. The matrix is symmetric for each group
(`pip[i,j,k] = pip[j,i,k]`). Diagonal entries are 0.

## See also

[`pip_diff_edge()`](https://mljaniczek.github.io/multiGGMr/reference/pip_diff_edge.md),
[`confusion_at_threshold()`](https://rdrr.io/pkg/ggmtools/man/confusion_at_threshold.html),
[`roc_auc()`](https://rdrr.io/pkg/ggmtools/man/roc_auc.html)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
pip <- pip_edges(fit)
# Edges with PIP >= 0.5 in group 1
sum(pip[,,1][upper.tri(pip[,,1])] >= 0.5)
#> [1] 8
```
