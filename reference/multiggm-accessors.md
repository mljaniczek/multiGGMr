# Shared-format accessors for multiGGMr fits

Methods for the ggmtools generics
[`ggmtools::pip_table()`](https://rdrr.io/pkg/ggmtools/man/pip_table.html)
and
[`ggmtools::pcor_estimate()`](https://rdrr.io/pkg/ggmtools/man/pcor_estimate.html),
so that downstream code written for those generics works on multiGGMr
fits as well as on fits from other packages (e.g. spikeyglass).

## Usage

``` r
# S3 method for class 'multiggm_fit'
pip_table(fit, varnames = NULL, ...)

# S3 method for class 'multiggm_fit_list'
pip_table(fit, chain = NULL, varnames = NULL, ...)

# S3 method for class 'multiggm_platform_fit'
pip_table(fit, platform, varnames = NULL, ...)

# S3 method for class 'multiggm_fit'
pcor_estimate(fit, varnames = NULL, ...)

# S3 method for class 'multiggm_fit_list'
pcor_estimate(fit, chain = NULL, varnames = NULL, ...)

# S3 method for class 'multiggm_platform_fit'
pcor_estimate(fit, platform, varnames = NULL, ...)
```

## Arguments

- fit:

  A `multiggm_fit`, `multiggm_fit_list` or `multiggm_platform_fit`.

- varnames:

  Optional variable names.

- ...:

  Ignored.

- chain:

  For a `multiggm_fit_list`: `NULL` (default) pools all chains; an
  integer selects chains.

- platform:

  Platform to use for a `multiggm_platform_fit`.

## Value

[`pip_table()`](https://rdrr.io/pkg/ggmtools/man/pip_table.html): a data
frame with columns `i`, `j`, `node1`, `node2`, `group`, `pip`.
[`pcor_estimate()`](https://rdrr.io/pkg/ggmtools/man/pcor_estimate.html):
a list of K partial-correlation matrices named `Group_1`, ...,
`Group_K`.

## Details

[`pip_table()`](https://rdrr.io/pkg/ggmtools/man/pip_table.html) returns
each group's own posterior edge-inclusion probabilities (from
[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md))
in long format. It needs edge indicators, so it stops for mGHS fits.
[`pcor_estimate()`](https://rdrr.io/pkg/ggmtools/man/pcor_estimate.html)
returns the posterior mean partial correlations (as
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)). For
multi-chain fits all chains are pooled unless `chain` selects one, as
for
[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md);
for multi-platform fits, choose `platform`.

multiGGMr fits do not store variable names, so the nodes are labelled
`V1, ..., Vp` unless `varnames` is given.

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 100, nsave = 50)
head(pip_table(fit))
#>   i j node1 node2   group  pip
#> 1 1 2    V1    V2 Group_1 0.46
#> 2 1 3    V1    V3 Group_1 0.98
#> 3 2 3    V2    V3 Group_1 1.00
#> 4 1 4    V1    V4 Group_1 0.00
#> 5 2 4    V2    V4 Group_1 0.66
#> 6 3 4    V3    V4 Group_1 0.30
pcor_estimate(fit)$Group_1[1:3, 1:3]
#>             V1          V2         V3
#> V1  1.00000000 -0.06305474 -0.2820263
#> V2 -0.06305474  1.00000000  0.3031916
#> V3 -0.28202632  0.30319164  1.0000000
```
