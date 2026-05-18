# Extract posterior summary of precision matrices

Computes a summary statistic (default: mean) of the precision matrix
draws for each group. More flexible than
[`coef.multiggm_fit`](https://mljaniczek.github.io/multiGGMr/reference/coef.multiggm_fit.md)
because you can specify an arbitrary summary function.

## Usage

``` r
posterior_precision(fit, summary_fun = mean, chain = 1L)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list` object.

- summary_fun:

  Function to apply across posterior draws for each matrix entry.
  Default `mean`. Other useful choices include `median` or `sd`.

- chain:

  Integer; which chain to use if `fit` is a `multiggm_fit_list`. Default
  1.

## Value

A named list of K numeric matrices (each p x p). List names are
`"Group_1"`, `"Group_2"`, etc.

## See also

[`coef.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/coef.multiggm_fit.md),
[`posterior_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
# Posterior mean
omega_mean <- posterior_precision(fit)
# Posterior standard deviation
omega_sd <- posterior_precision(fit, summary_fun = sd)
```
