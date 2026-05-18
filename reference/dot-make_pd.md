# Make a symmetric matrix positive definite

Rescales off-diagonal entries row-wise and adds diagonal loading if
needed.

## Usage

``` r
.make_pd(A, denom_factor = 1.5)
```

## Arguments

- A:

  Symmetric matrix with unit diagonal.

- denom_factor:

  Scale factor (larger = more shrinkage). Default 1.5.

## Value

A symmetric positive definite matrix.
