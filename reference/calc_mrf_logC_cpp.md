# Compute log normalizing constant for the edge-sharing MRF prior

Computes \\\log C(\Theta, \nu)\\ where
\$\$C(\Theta,\nu)=\sum\_{g\in\\0,1\\^K} \exp\\\nu \mathbf{1}^\top g +
g^\top \Theta g\\.\$\$

## Usage

``` r
calc_mrf_logC_cpp(Theta, nu)
```

## Arguments

- Theta:

  A K x K symmetric matrix.

- nu:

  Numeric vector of length 1 or more.

## Value

Numeric vector of \\\log C(\Theta,\nu)\\ for each entry of `nu`.
