# Extract posterior mean partial correlation matrices

Returns a list of K posterior mean partial correlation matrices from a
fitted multi-GGM model. This is the `fitted` S3 method for
`multiggm_fit` objects.

## Usage

``` r
# S3 method for class 'multiggm_fit'
fitted(object, ...)
```

## Arguments

- object:

  A `multiggm_fit` object returned by
  [`multiggm_mcmc`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md).

- ...:

  Ignored.

## Value

A named list of K numeric matrices (each p x p). Each matrix is the
element-wise posterior mean of the partial correlation matrix \\P_k\\
across all saved MCMC iterations, where \\P\_{ij}^{(k)} =
-\Omega\_{ij}^{(k)} / \sqrt{\Omega\_{ii}^{(k)} \Omega\_{jj}^{(k)}}\\.
Diagonal entries are 1. List names are `"Group_1"`, `"Group_2"`, etc.

## Details

For each saved iteration, the precision matrix \\\Omega_k\\ is converted
to a partial correlation matrix, then the posterior mean is taken
element-wise. This is computed via
[`posterior_pcor`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md).

## See also

[`coef.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/coef.multiggm_fit.md),
[`posterior_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md),
[`precision_to_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/precision_to_pcor.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
pcor_hat <- fitted(fit)
round(pcor_hat$Group_1[1:5, 1:5], 3)
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,]  1.000 -0.480 -0.014  0.005 -0.003
#> [2,] -0.480  1.000 -0.020  0.005  0.050
#> [3,] -0.014 -0.020  1.000  0.119 -0.097
#> [4,]  0.005  0.005  0.119  1.000 -0.219
#> [5,] -0.003  0.050 -0.097 -0.219  1.000
```
