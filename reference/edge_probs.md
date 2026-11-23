# Posterior edge probabilities for every group and pair of groups

Returns per-edge posterior summaries computed from running sums that
every sampler keeps over all post-burn-in iterations, so they are
available without saving draws (`save_draws = FALSE`). Each probability
is counted per MCMC iteration and then averaged, never assembled from
the groups' marginal summaries: the groups are dependent through the MRF
prior (G-Wishart, SSVS) or the correlation of edge scales (mGHS).

## Usage

``` r
edge_probs(fit, scale = c("precision", "pcor"), c = NULL, chain = NULL)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list`.

- scale:

  `"precision"` (default) or `"pcor"`: the scale of `mean`, `sd`,
  `p_exceed_*`, `p_union_*` and `p_sign`.

- c:

  Thresholds for the magnitude and union probabilities. `NULL` (default)
  uses the grid given by `diff_deltas` at fit time. Values outside that
  grid are computed from saved draws when the fit has them (for mGHS
  this also needs the shrinkage-weight draws), and are an error
  otherwise.

- chain:

  `NULL` (default) pools all chains; an integer selects chains.

## Value

A data frame with columns `type`, `group1`, `group2`, `i`, `j` and the
quantities above. Edges are in R's
[`upper.tri()`](https://rdrr.io/r/base/lower.tri.html) order.
Attributes: `scale`, `c`, `n_draws`, `method` and, for mGHS, `t_hat`.

## Details

Rows with `type = "group"` describe group `group1`; rows with
`type = "pair"` describe the difference \\D = x\_{group1} -
x\_{group2}\\ for every pair of groups. Columns:

- `mean`, `sd`:

  Posterior mean and SD of the precision entry or partial correlation
  (group rows) or of \\D\\ (pair rows).

- `p_struct`:

  Structural probability. Group rows: the edge's posterior inclusion
  probability \\P(G\_{k,ij} = 1)\\ for G-Wishart and SSVS; for mGHS
  \\P(w\_{k,ij} \> 0.5)\\, where \\w = \lambda^2 / (1 + \lambda^2)\\ is
  the shrinkage weight (the paper's \\\kappa\\). Pair rows: \\P(G_s \ne
  G_t)\\, or for mGHS \\P(S_s \ne S_t)\\ with \\S = 1(w \> 0.5)\\.

- `mean_weight`:

  mGHS group rows: \\E\[w\]\\. The authors' median-probability rule
  selects \\E\[w\] \> 0.5\\ (Busatto & Stingo, Section 5).

- `p_cut`:

  mGHS: probability of the cut-model indicator, \\P(Z = 1)\\ (group
  rows) or \\P(Z_s \ne Z_t)\\ (pair rows).

- `p_exceed_<c>`:

  \\P(\|x\| \> c)\\ (group rows) or \\P(\|D\| \> c)\\ (pair rows).

- `p_union_<c>`:

  Pair rows: \\P(\textrm{structural difference or } \|D\| \> c)\\,
  counted per iteration. The structural part is the edge indicator for
  G-Wishart and SSVS and the cut-model indicator for mGHS.

- `p_sign`:

  \\\max\\P(x \> 0), P(x \< 0)\\\\ (group rows) or \\\max\\P(D \> 0),
  P(D \< 0)\\\\ (pair rows). A central credible interval at level \\L\\
  excludes zero exactly when `p_sign` \\\> (1 + L) / 2\\; the graphical
  horseshoe's published rule (Li, Craig & Bhadra 2019) is \\L = 0.5\\,
  i.e. `p_sign > 0.75` (`select_edges(rule = "interval")`). For pairs it
  is **not** a probability that the edge differs (see
  [`diff_summary()`](https://mljaniczek.github.io/multiGGMr/reference/diff_summary.md)).

- `mean_shrinkage`, `p_shrinkage`, `mcse_shrinkage`:

  mGHS group rows: \\E\[1 - \kappa\]\\, \\P(1 - \kappa \> 0.5)\\ and the
  Monte Carlo SE of the latter, where \\\kappa\\ is the edge's shrinkage
  factor with the global scale included (see
  [`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md),
  `rule = "shrinkage"`).

- `p_up_<c>`, `p_down_<c>`:

  Pair rows: \\P(D \> c)\\ and \\P(D \< -c)\\, the direction of a
  difference: with \\D = x\_{group1} - x\_{group2}\\, `p_up` means the
  entry is larger in `group1`. They add up to `p_exceed_<c>`.

- `mcse_struct`, `mcse_cut`, `mcse_exceed_<c>`, `mcse_union_<c>`:

  Batch-means Monte Carlo standard errors of `p_struct`, `p_cut`,
  `p_exceed_<c>` and `p_union_<c>` (`NA` if the fit was made with
  `mcse_batches = 0` or by an older version). With several chains they
  are pooled as for the probabilities. An edge whose probability is
  within about two Monte Carlo SEs of a selection cutoff needs a longer
  run.

Columns that do not apply are `NA`.

**mGHS structural and cut-model probabilities are weak.** Following the
paper, \\w = \lambda^2 / (1 + \lambda^2)\\ uses the local scale only.
The horseshoe shrinks null edges mainly through the global scale
\\\tau\\, so for a null edge \\\lambda\\ stays close to its half-Cauchy
prior, under which \\P(w \> 0.5) = 0.5\\. Hence `p_struct` and
`mean_weight` are near 0.5 for null edges, and \\P(S_s \ne S_t)\\ is
near 0.5 for edges that are null in both groups. The cut-model indicator
is drawn as \\Z \sim \textrm{Bernoulli}(\Phi(w - t))\\ (the upstream
code; the paper uses the closed-form \\P(w \> t \mid \cdot)\\), so \\P(Z
= 1)\\ lies in about (0.16, 0.84) and \\P(Z_s \ne Z_t)\\ is near 0.5 for
most edges. For differential edges with mGHS, use the magnitude rule.
The authors select edges by \\E\[w\] \> 0.5\\
(`rule = "median_shrinkage"`) or by \\E\[w\] \> \hat t_k\\, with \\\hat
t_k\\ the posterior mean threshold (attribute `t_hat`).

For mGHS the result also carries the posterior mean cut-model threshold
per group as attribute `t_hat`, for the authors' rule \\E\[w\] \> \hat
t_k\\.

## See also

[`edge_prob_matrix()`](https://mljaniczek.github.io/multiGGMr/reference/edge_prob_matrix.md),
[`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md),
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
                     burnin = 100, nsave = 100, save_draws = FALSE)
ep <- edge_probs(fit)
head(ep[ep$type == "pair", ])
#>    type group1 group2 i j          mean         sd p_struct mean_weight p_cut
#> 31 pair      1      2 1 2  0.0141448044 0.03330287     0.05          NA    NA
#> 32 pair      1      2 1 3 -0.1462129269 0.21682931     0.42          NA    NA
#> 33 pair      1      2 2 3 -0.7640950740 0.20150045     0.08          NA    NA
#> 34 pair      1      2 1 4  0.0001849808 0.03113862     0.02          NA    NA
#> 35 pair      1      2 2 4 -0.1113987778 0.14902535     0.38          NA    NA
#> 36 pair      1      2 3 4 -0.0358687026 0.08926421     0.13          NA    NA
#>    p_sign p_exceed_0.01 p_exceed_0.02 p_exceed_0.03 p_exceed_0.05 p_exceed_0.1
#> 31   0.68          0.80          0.58          0.40          0.14         0.01
#> 32   0.74          0.98          0.95          0.95          0.91         0.78
#> 33   1.00          1.00          1.00          1.00          1.00         1.00
#> 34   0.50          0.72          0.45          0.26          0.08         0.01
#> 35   0.83          0.84          0.70          0.59          0.46         0.34
#> 36   0.67          0.74          0.50          0.38          0.20         0.11
#>    p_union_0.01 p_union_0.02 p_union_0.03 p_union_0.05 p_union_0.1 p_up_0.01
#> 31         0.80         0.58         0.41         0.16        0.05      0.57
#> 32         0.98         0.95         0.95         0.91        0.80      0.26
#> 33         1.00         1.00         1.00         1.00        1.00      0.00
#> 34         0.73         0.46         0.27         0.09        0.02      0.36
#> 35         0.86         0.72         0.61         0.49        0.42      0.10
#> 36         0.75         0.51         0.39         0.21        0.15      0.23
#>    p_up_0.02 p_up_0.03 p_up_0.05 p_up_0.1 p_down_0.01 p_down_0.02 p_down_0.03
#> 31      0.42      0.32      0.12     0.01        0.23        0.16        0.08
#> 32      0.25      0.25      0.22     0.17        0.72        0.70        0.70
#> 33      0.00      0.00      0.00     0.00        1.00        1.00        1.00
#> 34      0.21      0.12      0.03     0.01        0.36        0.24        0.14
#> 35      0.06      0.05      0.00     0.00        0.74        0.64        0.54
#> 36      0.13      0.09      0.03     0.00        0.51        0.37        0.29
#>    p_down_0.05 p_down_0.1 mcse_struct mcse_cut mcse_exceed_0.01
#> 31        0.02       0.00  0.03516877       NA       0.03554093
#> 32        0.69       0.61  0.10042017       NA       0.01376494
#> 33        1.00       1.00  0.05109331       NA       0.00000000
#> 34        0.05       0.00  0.02000000       NA       0.05505978
#> 35        0.46       0.34  0.09049571       NA       0.04000000
#> 36        0.17       0.11  0.05482028       NA       0.03276712
#>    mcse_exceed_0.02 mcse_exceed_0.03 mcse_exceed_0.05 mcse_exceed_0.1
#> 31       0.05005260       0.05428967       0.03866183      0.01000000
#> 32       0.02460210       0.02460210       0.03692097      0.05211324
#> 33       0.00000000       0.00000000       0.00000000      0.00000000
#> 34       0.05000000       0.04834089       0.02247806      0.01000000
#> 35       0.06407233       0.06878647       0.08220642      0.08472867
#> 36       0.04920419       0.06307890       0.06488857      0.05124759
#>    mcse_union_0.01 mcse_union_0.02 mcse_union_0.03 mcse_union_0.05
#> 31      0.03554093      0.05005260      0.05520297      0.04723959
#> 32      0.01376494      0.02460210      0.02460210      0.03692097
#> 33      0.00000000      0.00000000      0.00000000      0.00000000
#> 34      0.05482028      0.05047146      0.05083513      0.02704772
#> 35      0.03866183      0.06390783      0.07030011      0.08642003
#> 36      0.03203616      0.04915068      0.06236565      0.06403124
#>    mcse_union_0.1 mean_shrinkage p_shrinkage mcse_shrinkage
#> 31     0.03516877             NA          NA             NA
#> 32     0.05026247             NA          NA             NA
#> 33     0.00000000             NA          NA             NA
#> 34     0.02000000             NA          NA             NA
#> 35     0.09049571             NA          NA             NA
#> 36     0.06303717             NA          NA             NA
```
