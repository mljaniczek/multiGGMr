# Rescale a symmetric matrix toward positive definiteness

Port of the helper in the MATLAB reference code. For each row, rescales
off-diagonal entries by `(denom_factor * row L1 norm of off-diagonals)`,
then enforces unit diagonal and symmetry. This makes the matrix
diagonally dominant (hence positive definite) when `denom_factor > 1`.

## Usage

``` r
fix_matrix(A, denom_factor = 1)
```

## Arguments

- A:

  Numeric square matrix with unit diagonal.

- denom_factor:

  Positive scalar controlling shrinkage of off-diagonal entries. Values
  \> 1 shrink off-diagonals more aggressively. Default 1.

## Value

A symmetric numeric matrix with unit diagonal, where each row's
off-diagonal entries have been rescaled so their absolute sum is at most
`1 / denom_factor`.

## Examples

``` r
A <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
fix_matrix(A, denom_factor = 1.5)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.6666667
#> [2,] 0.6666667 1.0000000
```
