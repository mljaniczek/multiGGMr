# Summarize a multi-chain multiggm fit

Summarizes each chain and the pooled posterior. The pooled edge counts
use posterior inclusion probabilities averaged over all chains. The
largest difference between two chains' PIPs (or, for mGHS, posterior
mean partial correlations) is a quick convergence check: chains that
agree give small values.

## Usage

``` r
# S3 method for class 'multiggm_fit_list'
summary(object, pip_threshold = 0.5, ...)
```

## Arguments

- object:

  A `multiggm_fit_list` returned by
  [`multiggm_mcmc`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md)
  with `nchains > 1`.

- pip_threshold:

  Numeric; PIP threshold for the edge counts. Default 0.5.

- ...:

  Ignored.

## Value

An object of class `"summary.multiggm_fit_list"`: a list with `chains`
(the per-chain summaries), `edge_counts` (pooled, per group, or `NULL`
for mGHS), `max_chain_diff` (largest absolute between-chain difference,
per group), `nchains`, `K` and `pip_threshold`.

## See also

[`summary.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/summary.multiggm_fit.md),
[`plot.multiggm_fit_list()`](https://mljaniczek.github.io/multiGGMr/reference/plot.multiggm_fit.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 100, nsave = 50,
                     nchains = 2, seed = 1)
summary(fit)
#> multiGGM MCMC Summary: 2 chains, K = 2
#> ==========================================
#> Pooled selected edges (PIP >= 0.5 ):
#>   Group 1: 7 edges
#>   Group 2: 7 edges
#> 
#> Largest difference in PIP between chains:
#>   Group 1: 0.480
#>   Group 2: 0.520
#> 
#> --- Chain 1 ---
#> multiGGM MCMC Summary
#> =====================
#> Groups (K): 2   |  Nodes (p): 6   |  Posterior draws: 50 
#> 
#> Acceptance rates:
#>   gamma (edge toggle): 6.0% 
#>   theta (within-model): 36.3% 
#>   nu (edge log-odds): 58.5% 
#> 
#> Selected edges (PIP >= 0.5 ):
#>   Group 1 : 8 edges
#>   Group 2 : 8 edges
#> 
#> Graph similarity (theta):
#>   theta[1,2]: mean = 1.022, P(nonzero) = 98.0%
#> 
#> 
#> --- Chain 2 ---
#> multiGGM MCMC Summary
#> =====================
#> Groups (K): 2   |  Nodes (p): 6   |  Posterior draws: 50 
#> 
#> Acceptance rates:
#>   gamma (edge toggle): 2.0% 
#>   theta (within-model): 39.6% 
#>   nu (edge log-odds): 59.4% 
#> 
#> Selected edges (PIP >= 0.5 ):
#>   Group 1 : 7 edges
#>   Group 2 : 7 edges
#> 
#> Graph similarity (theta):
#>   theta[1,2]: mean = 0.764, P(nonzero) = 98.0%
#> 
```
