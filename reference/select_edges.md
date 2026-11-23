# Select edges within one group by posterior probability

Selects the edges of one group's graph from the per-edge posterior
probabilities of
[`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md),
either by the median-probability rule (probability \> 0.5) or by
Bayesian FDR control
([`bayes_fdr()`](https://mljaniczek.github.io/multiGGMr/reference/bayes_fdr.md)).

## Usage

``` r
select_edges(
  fit,
  group = 1,
  rule = c("structural", "magnitude", "cut", "median_shrinkage", "shrinkage", "interval"),
  c = NULL,
  decision = NULL,
  alpha = 0.05,
  scale = c("precision", "pcor"),
  chain = NULL,
  level = 0.5
)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list`.

- group:

  The group whose graph to select. Default 1.

- rule:

  `"structural"`

  :   (default) \\P(G\_{k,ij} = 1)\\ (G-Wishart, SSVS) or \\P(w\_{k,ij}
      \> 0.5)\\ (mGHS).

  `"magnitude"`

  :   \\P(\|x\_{k,ij}\| \> c)\\ on `scale`.

  `"cut"`

  :   mGHS only: \\P(Z\_{k,ij} = 1)\\, the cut-model indicator.

  `"median_shrinkage"`

  :   mGHS only: the authors' median rule \\E\[w\_{k,ij}\] \> 0.5\\.
      \\E\[w\]\\ is not a probability, so only `decision = "median"`
      applies and no FDR is estimated.

  `"shrinkage"`

  :   mGHS only: \\E\[1 - \kappa\_{k,ij}\]\\, where \\\kappa\\ is the
      edge's shrinkage factor in its column update with the global scale
      included: with \\a = \mu \tau_k \lambda^2\_{ij,k} s\_{jj}
      (\Omega\_{-j}^{-1})\_{ii}\\ (prior variance times likelihood
      precision), \\1 - \kappa = a / (1 + a)\\, averaged over the two
      column updates of the edge and over iterations. In horseshoe
      regression \\1 - \kappa\\ is read as a pseudo inclusion
      probability (Carvalho, Polson & Scott 2010). It measures how much
      the data dominate the prior for the edge, not whether the edge is
      nonzero, so it grows with \\n\\ for null edges too. With
      correlated groups an edge is shrunk towards the value the other
      groups predict, and the approximation ignores coupling between
      edges. `decision = "median"` (\\E\[1 - \kappa\] \> 0.5\\) is the
      default; `"bfdr"` treats it as a probability. In simulations (p =
      8 and 20, n = 100 to 1000) it ranked edges as well as \\E\[w\]\\,
      but neither 0.5 rule controlled false edges: at p = 20, n = 100,
      \\E\[1 - \kappa\] \> 0.5\\ selected 3 null edges against 65 for
      \\E\[w\] \> 0.5\\; at p = 8, n = 300, 17 against 0. For selection
      with error control use `rule = "magnitude"`.

  `"interval"`

  :   The central credible interval at level `level` excludes zero, i.e.
      `p_sign` \\\> (1 + \\ `level`\\) / 2\\. With `level = 0.5` this is
      the graphical horseshoe's published rule (Li, Craig & Bhadra
      2019). It is a sign rule: as \\n\\ grows it selects null edges
      whose small estimates are well determined, so it is offered for
      comparison, with no FDR estimate.

- c:

  Threshold for `rule = "magnitude"`.

- decision:

  `"bfdr"` (default; `"median"` for `rule = "median_shrinkage"` and
  `"shrinkage"`) selects the largest set with posterior expected FDR at
  most `alpha`, keeping tied probabilities together; `"median"` selects
  probability \> 0.5.

- alpha:

  Target posterior expected FDR for `decision = "bfdr"`.

- scale:

  `"precision"` (default) or `"pcor"`: the scale of `mean`, `sd`,
  `p_exceed_*`, `p_union_*` and `p_sign`.

- chain:

  `NULL` (default) pools all chains; an integer selects chains.

- level:

  Credible level for `rule = "interval"` (default 0.5).

## Value

An object of class `"multiggm_edge_selection"`: a list with `selected`
(logical `p x p`), `prob`, `edges` (the selected edges with their
probability and Bayesian FDR q-value, see
[`bfdr_qvalues()`](https://mljaniczek.github.io/multiGGMr/reference/bfdr_qvalues.md)),
`expected_fdr` (posterior expected FDR of the selected set),
`n_selected`, `alpha`, `rule`, `c`, `decision`, `scale` and `group`.

## See also

[`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md),
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md),
[`bayes_fdr()`](https://mljaniczek.github.io/multiGGMr/reference/bayes_fdr.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
                     burnin = 100, nsave = 100, save_draws = FALSE)
select_edges(fit, group = 1, decision = "median")
#> <multiggm_edge_selection>
#>   Group 1, rule: structural, decision: median
#>   Selected 4 edges with value > 0.5 (posterior expected FDR 0.175)
#>   i j prob    q_value
#> 1 5 6 1.00 0.00000000
#> 2 2 3 0.92 0.04000000
#> 3 4 5 0.80 0.09333333
#> 4 1 3 0.58 0.17500000
```
