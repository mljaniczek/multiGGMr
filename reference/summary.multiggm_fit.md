# Summarize a multiggm_fit object

Displays MCMC diagnostics including acceptance rates, edge counts at a
given PIP threshold, and graph similarity (theta) estimates.

## Usage

``` r
# S3 method for class 'multiggm_fit'
summary(object, pip_threshold = 0.5, ...)
```

## Arguments

- object:

  A `multiggm_fit` object returned by
  [`multiggm_mcmc`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md).

- pip_threshold:

  Numeric threshold for counting selected edges via posterior inclusion
  probability (PIP). Default 0.5.

- ...:

  Ignored.

## Value

An object of class `"summary.multiggm_fit"` (printed invisibly) with
components:

- `K`:

  Integer; number of groups.

- `p`:

  Integer; number of variables.

- `nsave`:

  Integer; number of saved posterior draws.

- `ar_gamma`:

  Numeric; mean acceptance rate for between-model (spike-slab toggle)
  moves on theta.

- `ar_theta`:

  Numeric; mean acceptance rate for within-model (slab-to-slab) moves on
  theta.

- `ar_nu`:

  Numeric; mean acceptance rate for nu updates.

- `edge_counts`:

  Integer vector of length K; number of selected edges per group at the
  given PIP threshold.

- `pip_threshold`:

  The PIP threshold used.

- `theta_mean`:

  Numeric matrix `[K, K]`; posterior mean of theta (graph similarity)
  for each pair of groups. Upper triangle only.

- `theta_nonzero_frac`:

  Numeric matrix `[K, K]`; posterior probability that theta \> 0 for
  each pair. Upper triangle only.

- `hyper`:

  Named list of hyperparameters used in the fit.

## See also

[`multiggm_mcmc()`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md),
[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
summary(fit)
#> multiGGM MCMC Summary
#> =====================
#> Groups (K): 2   |  Nodes (p): 8   |  Posterior draws: 100 
#> 
#> Acceptance rates:
#>   gamma (edge toggle): 7.7% 
#>   theta (within-model): 62.7% 
#>   nu (edge log-odds): 56.0% 
#> 
#> Selected edges (PIP >= 0.5 ):
#>   Group 1 : 8 edges
#>   Group 2 : 7 edges
#> 
#> Graph similarity (theta):
#>   theta[1,2]: mean = 0.493, P(nonzero) = 97.0%
#> 
summary(fit, pip_threshold = 0.3)
#> multiGGM MCMC Summary
#> =====================
#> Groups (K): 2   |  Nodes (p): 8   |  Posterior draws: 100 
#> 
#> Acceptance rates:
#>   gamma (edge toggle): 7.7% 
#>   theta (within-model): 62.7% 
#>   nu (edge log-odds): 56.0% 
#> 
#> Selected edges (PIP >= 0.3 ):
#>   Group 1 : 11 edges
#>   Group 2 : 10 edges
#> 
#> Graph similarity (theta):
#>   theta[1,2]: mean = 0.493, P(nonzero) = 97.0%
#> 
```
