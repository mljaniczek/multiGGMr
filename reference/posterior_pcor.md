# Compute posterior partial correlations from a multiggm_fit object

Extracts the full array of posterior partial correlation draws from a
fitted multi-GGM model. For each saved MCMC iteration and group, the
precision matrix is converted to a partial correlation matrix.

## Usage

``` r
posterior_pcor(fit, use = c("C_save", "Sig_save"), chain = 1L)
```

## Arguments

- fit:

  A `multiggm_fit` object (single chain) or `multiggm_fit_list` object
  (multiple chains).

- use:

  Character; which posterior array to compute partial correlations from:

  `"C_save"`

  :   (default) Use the precision matrix draws directly. Faster and
      recommended.

  `"Sig_save"`

  :   Use the covariance matrix draws, inverting each to get precision
      first. Slower; mainly for cross-checking.

- chain:

  Integer; which chain to use if `fit` is a `multiggm_fit_list`. Default
  1.

## Value

A numeric array of dimension `[p, p, K, nsave]`. Entry
`pcor[i, j, k, s]` is the partial correlation between variables `i` and
`j` in group `k` at saved iteration `s`. Diagonal entries are 1.
Off-diagonal entries are in `[-1, 1]`.

## Details

The partial correlation between variables `i` and `j` given all other
variables is defined as: \$\$\rho\_{ij\|rest} = -\Omega\_{ij} /
\sqrt{\Omega\_{ii} \Omega\_{jj}}\$\$ where \\\Omega\\ is the precision
matrix.

## See also

[`precision_to_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/precision_to_pcor.md),
[`posterior_ci()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_ci.md),
[`fitted.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/fitted.multiggm_fit.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
pcor_draws <- posterior_pcor(fit)
# Posterior mean partial correlation for group 1
pcor_mean <- apply(pcor_draws[,,1,], c(1,2), mean)
```
