# WAIC from running sums

The widely applicable information criterion (Watanabe 2010) of the
Gaussian likelihood, computed from per-observation running sums kept
during sampling, so it needs no saved draws. Fit with `waic = TRUE` and
`data_list`.

## Usage

``` r
compute_waic(fit, chain = NULL)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list` fitted with `waic = TRUE`.

- chain:

  `NULL` (default) pools all chains.

## Value

An object of class `"multiggm_waic"`: a list with `estimates` (data
frame with rows `elpd_waic`, `p_waic`, `waic` and columns `estimate`,
`se`), `by_group` (the same per group) and `pointwise` (data frame with
`group`, `obs`, `elpd_waic`, `p_waic`).

## Details

For observation \\i\\ with log-likelihood \\\ell\_{is}\\ at iteration
\\s = 1, \dots, S\\: \\\textrm{lppd}\_i = \log(S^{-1} \sum_s
e^{\ell\_{is}})\\, \\p_i = \textrm{var}\_s(\ell\_{is})\\ and
\\\textrm{elpd}\_i = \textrm{lppd}\_i - p_i\\. WAIC is \\-2 \sum_i
\textrm{elpd}\_i\\, lower is better. Use it to compare fits to the same
`data_list`, for example different `hyper` settings, samplers, or a
joint fit against separate fits per group (compare the per-group rows).
The data are column-centered as for fitting. Observations with \\p_i \>
0.4\\ make WAIC unreliable (a warning is given); PSIS-LOO from the loo
package is then preferable and needs saved draws.

## References

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
*Journal of Machine Learning Research*, 11, 3571-3594.

Vehtari, A., Gelman, A. & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27, 1413-1432.

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
                     burnin = 100, nsave = 100, save_draws = FALSE, waic = TRUE)
compute_waic(fit)
#> Warning: 21 of 120 observations have p_waic > 0.4; WAIC may be unreliable.
#> <multiggm_waic> computed from 120 observations
#>           estimate   se
#> elpd_waic  -1054.5 18.5
#> p_waic        28.3  2.6
#> waic        2109.0 37.0
```
