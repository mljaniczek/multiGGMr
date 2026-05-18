# Convert precision matrix to partial correlation matrix

Given a precision matrix \\\Omega\\, returns the partial correlation
matrix \\P\\ with entries \\P\_{ij} =
-\Omega\_{ij}/\sqrt{\Omega\_{ii}\\\Omega\_{jj}}\\ and ones on the
diagonal. Partial correlations measure the linear association between
variables `i` and `j` after conditioning on all other variables.

## Usage

``` r
precision_to_pcor(Omega)
```

## Arguments

- Omega:

  A square numeric matrix (precision / inverse covariance). Must have
  positive diagonal entries.

## Value

A symmetric numeric matrix of the same dimension as `Omega`, with
diagonal entries equal to 1 and off-diagonal entries in `[-1, 1]`.

## See also

[`posterior_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md),
[`fitted.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/fitted.multiggm_fit.md)

## Examples

``` r
Omega <- matrix(c(2, -0.5, -0.5, 1.5), 2, 2)
precision_to_pcor(Omega)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.2886751
#> [2,] 0.2886751 1.0000000
```
