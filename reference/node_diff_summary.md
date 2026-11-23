# Node-level differential network summaries

For each node, posterior summaries of how much its connections differ
between two groups, counted per MCMC iteration: the expected number of
differential edges at the node (its differential degree), the
probability that it has at least one, and its expected rewiring score
\\\sum_j \|D\_{ij}\|\\. These identify hub changes that per-edge
selection can miss when many small changes meet at one node.

## Usage

``` r
node_diff_summary(
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

  Threshold for `"magnitude"` and `"union"`; must be in the fit's
  `diff_deltas` grid unless draws were saved.

- scale:

  Scale of \\D\\ and of the rewiring score: `"pcor"` (default) or
  `"precision"`.

- chain:

  `NULL` (default) pools all chains; an integer vector selects chains of
  a `multiggm_fit_list`.

## Value

A data frame with one row per node: `node`, `diff_degree` (posterior
mean number of differential edges), `p_any` (posterior probability of at
least one), `rewiring` (posterior mean of \\\sum_j \|D\_{ij}\|\\ on
`scale`) and `rewiring_rank`. Attributes `groups`, `rule`, `c` and
`scale`.

## See also

[`diff_count_posterior()`](https://mljaniczek.github.io/multiGGMr/reference/diff_count_posterior.md),
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
                     burnin = 100, nsave = 100, save_draws = FALSE)
node_diff_summary(fit, c = 0.1)
#>   node diff_degree p_any  rewiring rewiring_rank
#> 1    1        0.71  0.69 0.3103390             6
#> 2    2        1.62  1.00 0.9619459             3
#> 3    3        2.80  1.00 1.3453996             1
#> 4    4        1.91  0.98 0.7001043             5
#> 5    5        2.78  1.00 1.2980500             2
#> 6    6        1.96  1.00 0.8990888             4
```
