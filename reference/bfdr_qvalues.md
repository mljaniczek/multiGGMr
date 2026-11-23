# Bayesian FDR q-values

For each edge, the smallest target `alpha` at which
[`bayes_fdr()`](https://mljaniczek.github.io/multiGGMr/reference/bayes_fdr.md)
would select it: the posterior expected FDR of the set of edges with
probability at least as high (ties kept together). Selecting the edges
with `q <= alpha` gives exactly the selection of
`bayes_fdr(prob, alpha)`.

## Usage

``` r
bfdr_qvalues(prob)
```

## Arguments

- prob:

  A vector of posterior probabilities, or a square matrix (the upper
  triangle is used).

## Value

Values of the same shape as `prob` (a symmetric matrix with `NA`
diagonal for matrix input).

## See also

[`bayes_fdr()`](https://mljaniczek.github.io/multiGGMr/reference/bayes_fdr.md)

## Examples

``` r
bfdr_qvalues(c(0.99, 0.95, 0.6, 0.2))
#> [1] 0.0100000 0.0300000 0.1533333 0.3150000
```
