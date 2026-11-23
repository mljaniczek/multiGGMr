# Posterior distribution of the number of differential edges

The number of edges that differ between two groups, recorded at every
post-burn-in iteration. Its posterior gives a network-level answer to
"do these networks differ?" (`p_any`, the posterior probability of at
least one differential edge) and an estimate with uncertainty of how
much they differ, without choosing an FDR level.

## Usage

``` r
diff_count_posterior(
  fit,
  groups = c(1, 2),
  rule = c("magnitude", "union", "structural", "cut"),
  c = NULL,
  scale = c("pcor", "precision"),
  chain = NULL
)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list` object.

- groups:

  Integer vector of length 2, the groups `c(k, m)` to compare. Default
  `c(1, 2)`.

- rule:

  Which edges count as differential: `"magnitude"` (\\\|D\| \> c\\,
  default), `"union"` (structural difference or \\\|D\| \> c\\),
  `"structural"` or `"cut"` (mGHS).

- c:

  Threshold for `"magnitude"` and `"union"`.

- scale:

  Scale of \\D\\ and of the rewiring score: `"pcor"` (default) or
  `"precision"`.

- chain:

  `NULL` (default) pools all chains; an integer vector selects chains of
  a `multiggm_fit_list`.

## Value

An object of class `"multiggm_diff_count"`: a list with `draws` (the
count at every iteration, chains concatenated), `mean`, `sd`,
`quantiles` (2.5%, 50%, 97.5%), `p_any`, `groups`, `rule`, `c` and
`scale`.

## Details

With continuous priors (SSVS spike, mGHS) every \\D\\ is nonzero, so the
count is meaningful only with a threshold `c` above sampling noise (see
the edge-selection vignette). Compare the posterior with the same count
under the prior (fit with `prior_only = TRUE`). The posterior odds of
`p_any` divided by its prior odds is a network-level Bayes factor for
any difference, informative only when the prior probability is not close
to 1 (with many edges and a small `c` the prior almost surely has some
difference).

## See also

[`node_diff_summary()`](https://mljaniczek.github.io/multiGGMr/reference/node_diff_summary.md),
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
                     burnin = 100, nsave = 100, save_draws = FALSE)
diff_count_posterior(fit, rule = "union", c = 0.1)
#> <multiggm_diff_count> groups 1 vs 2 | rule: union (c = 0.1, pcor scale)
#>   Differential edges: mean 6.2, sd 0.9, 95% interval [5, 8]
#>   P(at least one differential edge) = 1.000
```
