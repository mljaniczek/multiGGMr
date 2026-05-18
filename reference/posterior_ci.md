# Posterior credible intervals for partial correlations

Computes pointwise posterior credible intervals for each partial
correlation entry across groups.

## Usage

``` r
posterior_ci(pcor, probs = c(0.025, 0.975))
```

## Arguments

- pcor:

  A numeric array of dimension `[p, p, K, nsave]`, as returned by
  [`posterior_pcor`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md).

- probs:

  Numeric vector of length 2 giving the lower and upper quantile
  probabilities. Default `c(0.025, 0.975)` for 95\\ credible intervals.

## Value

A list with four components:

- `lower`:

  Numeric array `[p, p, K]`; lower quantile of each partial correlation.

- `median`:

  Numeric array `[p, p, K]`; posterior median of each partial
  correlation.

- `upper`:

  Numeric array `[p, p, K]`; upper quantile of each partial correlation.

- `probs`:

  The quantile probabilities used.

## See also

[`posterior_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md),
[`diff_prob_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/diff_prob_pcor.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
pcor_draws <- posterior_pcor(fit)
ci <- posterior_ci(pcor_draws)
# 95% CI for edge (1,2) in group 1
cat(sprintf("pcor[1,2]: %.3f (%.3f, %.3f)\n",
    ci$median[1,2,1], ci$lower[1,2,1], ci$upper[1,2,1]))
#> pcor[1,2]: -0.489 (-0.634, -0.327)
```
