# Log normalizing constant for the MRF prior (per edge)

Computes the log normalizing constant \\\log C(\Theta, \nu)\\ for the
Markov random field prior on graph edge indicators, as defined in
Peterson et al. (2015), Equation 3.2: \$\$C(\Theta, \nu) = \sum\_{g \in
\\0,1\\^K} \exp(\nu \cdot \mathbf{1}'g + g'\Theta g)\$\$ This is the
computational hot spot in the MRF prior evaluation and is implemented in
C++ via RcppArmadillo for performance.

## Usage

``` r
calc_mrf_logC(Theta, nu)
```

## Arguments

- Theta:

  Numeric matrix (K x K); symmetric graph similarity parameters.
  Off-diagonal entry \\\Theta\_{km}\\ controls the tendency for edge
  inclusion to be shared between groups k and m.

- nu:

  Numeric vector of length 1 or more; edge-specific log-odds parameters
  \\\nu\_{ij}\\. Each value controls the baseline edge inclusion
  probability for one edge across all groups.

## Value

Numeric vector of the same length as `nu`, where each entry is \\\log
C(\Theta, \nu\_{ij})\\ for the corresponding \\\nu\_{ij}\\.

## Details

The computation requires summing over all \\2^K\\ binary vectors in
\\\\0,1\\^K\\, which is feasible for moderate K (the method is designed
for K up to about 20). For each binary vector \\g\\, the term \\\nu
\cdot \mathbf{1}'g + g'\Theta g\\ is computed and the log-sum-exp trick
is used for numerical stability.

## See also

[`multiggm_mcmc()`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md)

## Examples

``` r
# Two groups, no similarity
Theta <- matrix(0, 2, 2)
calc_mrf_logC(Theta, nu = c(-1, 0, 1))
#> [1] 0.6265234 1.3862944 2.6265234

# Two groups with similarity theta = 0.5
Theta[1,2] <- Theta[2,1] <- 0.5
calc_mrf_logC(Theta, nu = c(-1, 0, 1))
#> [1] 0.7436684 1.7436684 3.2779784
```
