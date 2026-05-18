# Extract posterior mean precision matrices

Returns a list of K posterior mean precision matrices from a fitted
multi-GGM model. This is the `coef` S3 method for `multiggm_fit`
objects.

## Usage

``` r
# S3 method for class 'multiggm_fit'
coef(object, ...)
```

## Arguments

- object:

  A `multiggm_fit` object returned by
  [`multiggm_mcmc`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md).

- ...:

  Ignored.

## Value

A named list of K numeric matrices (each p x p). Each matrix is the
element-wise posterior mean of the precision matrix \\\Omega_k\\ across
all saved MCMC iterations. List names are `"Group_1"`, `"Group_2"`, etc.

## Details

The precision matrix \\\Omega_k\\ encodes the conditional independence
structure of group `k`: \\\Omega\_{ij}^{(k)} = 0\\ if and only if
variables `i` and `j` are conditionally independent in group `k`. The
posterior mean is computed from `object$C_save`.

For partial correlations (standardized precision), use
[`fitted.multiggm_fit`](https://mljaniczek.github.io/multiGGMr/reference/fitted.multiggm_fit.md)
instead.

## See also

[`fitted.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/fitted.multiggm_fit.md),
[`posterior_precision()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_precision.md),
[`posterior_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
omega_hat <- coef(fit)
round(omega_hat$Group_1[1:5, 1:5], 3)
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,]  1.605  0.642  0.016 -0.008  0.003
#> [2,]  0.642  1.102  0.021 -0.007 -0.054
#> [3,]  0.016  0.021  0.858 -0.132  0.096
#> [4,] -0.008 -0.007 -0.132  1.390  0.267
#> [5,]  0.003 -0.054  0.096  0.267  1.031
```
