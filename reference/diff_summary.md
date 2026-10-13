# Per-edge posterior summary of a group difference

Summarizes the posterior of \\D\_{ij} = X\_{k,ij} - X\_{m,ij}\\ for
every edge \\i \< j\\: mean, standard deviation, credible interval, the
sign probability and the exceedance probabilities \\P(\|D\_{ij}\| \>
c)\\.

## Usage

``` r
diff_summary(
  fit,
  groups = c(1, 2),
  scale = c("pcor", "precision"),
  deltas = c(0.01, 0.02, 0.05, 0.1),
  probs = c(0.025, 0.975),
  chain = NULL
)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list` object.

- groups:

  Integer vector of length 2, the groups `c(k, m)` to compare. Default
  `c(1, 2)`.

- scale:

  `"pcor"` (default) for partial correlations or `"precision"` for
  precision-matrix entries.

- deltas:

  Numeric vector of magnitude thresholds `c` for the exceedance
  probabilities \\P(\|D\| \> c)\\, on the scale given by `scale`.

- probs:

  Numeric vector of length 2 with the credible-interval quantiles, or
  `NULL` to skip the interval (faster).

- chain:

  `NULL` (default) pools all chains; an integer vector selects chains of
  a `multiggm_fit_list`.

## Value

A data frame with one row per edge \\i \< j\\ and columns `i`, `j`,
`mean`, `sd`, `lower`, `median`, `upper`, `excludes_zero`, `p_sign`, and
one `p_exceed_<c>` column per entry of `deltas`. Rows are in R's
[`upper.tri()`](https://rdrr.io/r/base/lower.tri.html) order. Attributes
record `groups`, `scale`, `deltas`, `probs` and `n_draws`.

## Details

**Do not use `excludes_zero` or `p_sign` for FDR control.** A credible
interval that excludes zero is equivalent to `p_sign` \\\ge 1 -
\alpha/2\\, where `p_sign` \\= \max\\P(D \> 0), P(D \< 0)\\\\. Whenever
the posterior of \\D\\ is continuous, \\P(D = 0) = 0\\, so `p_sign` is
not a probability that the edge differs. This holds for the SSVS sampler
(normal spike), for the mGHS sampler (continuous shrinkage), and for
G-Wishart edges present in both groups. Edges with identical true values
in both groups then get `p_sign` near 1 from sampling noise alone. In a
simulation with mGHS (p = 100, n = 1000 per group), Bayesian FDR on
`p_sign` gave an achieved FDR of about 0.50 at a nominal 0.05. Use
\\P(\|D\| \> c)\\ with a practically meaningful `c` instead (see
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)).

For mGHS fits with `save_draws = FALSE`, the summary is read from
running sums kept by the sampler over every post-burn-in iteration. Then
`deltas` must be a subset of the grid given at fit time, and the
credible-interval columns are `NA`.

## See also

[`diff_draws()`](https://mljaniczek.github.io/multiGGMr/reference/diff_draws.md),
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md),
[`bayes_fdr()`](https://mljaniczek.github.io/multiGGMr/reference/bayes_fdr.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
ds <- diff_summary(fit, deltas = c(0.05, 0.1))
head(ds)
#>   i j         mean         sd        lower     median      upper excludes_zero
#> 1 1 2 -0.545838605 0.11333604 -0.762861810 -0.5249929 -0.3813812          TRUE
#> 2 1 3  0.007934572 0.08006625 -0.178854518  0.0000000  0.1858638         FALSE
#> 3 2 3  0.019830943 0.10398023 -0.154502805  0.0000000  0.2663474         FALSE
#> 4 1 4  0.008392163 0.04324758 -0.007186162  0.0000000  0.1475525         FALSE
#> 5 2 4  0.419292562 0.09993498  0.245943073  0.4058617  0.6534092          TRUE
#> 6 3 4  0.499296070 0.16104462  0.239851834  0.4925896  0.7949466          TRUE
#>   p_sign p_exceed_0.05 p_exceed_0.1
#> 1   1.00          1.00         1.00
#> 2   0.15          0.16         0.14
#> 3   0.21          0.27         0.19
#> 4   0.09          0.08         0.04
#> 5   1.00          1.00         1.00
#> 6   0.98          0.98         0.98
```
