# Prior versus posterior edge probabilities and Bayes factors

Compares each edge probability from
[`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md)
with the same probability under the prior alone, estimated from a fit
with `prior_only = TRUE`. The prior probability of an edge, or of a
difference \\\|D\| \> c\\, is set by the MRF and slab/spike or G-Wishart
priors and is generally not 0.5, so a posterior probability is only
interpretable next to it. The Bayes factor (posterior odds / prior odds)
measures how much the data moved each edge.

## Usage

``` r
edge_bayes_factors(
  fit,
  prior_fit,
  quantity = c("p_struct", "p_exceed", "p_union"),
  c = NULL,
  scale = c("pcor", "precision"),
  pool_edges = TRUE
)
```

## Arguments

- fit:

  The fitted model.

- prior_fit:

  A fit of the same model (same `method`, `hyper`, groups, `p` and
  `diff_deltas`) with `prior_only = TRUE`.

- quantity:

  `"p_struct"` (default), `"p_exceed"` or `"p_union"`.

- c:

  Threshold for `"p_exceed"` and `"p_union"`.

- scale:

  Scale of \\x\\ or \\D\\ for the magnitude quantities.

- pool_edges:

  If `TRUE` (default), the prior probability is averaged over the edges
  of each group (or pair). The default priors treat all edges alike, so
  this removes Monte Carlo noise; use `FALSE` for edge-specific priors.

## Value

A data frame with `type`, `group1`, `group2`, `i`, `j`, `prior`,
`posterior`, `bayes_factor` and `log10_bf`.

## Details

To keep Bayes factors finite, probabilities are estimated as \\(x + 0.5)
/ (n + 1)\\ from \\x\\ of \\n\\ iterations, which bounds the Bayes
factor by about \\2n\\ times the prior odds. Not available for mGHS,
whose diagonal prior is improper; under its half-Cauchy local scales the
shrinkage weight has a Beta(1/2, 1/2) prior and \\P(w \> 0.5) = 0.5\\.

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
                     burnin = 100, nsave = 200, save_draws = FALSE, seed = 1)
prior <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
                       burnin = 100, nsave = 200, save_draws = FALSE,
                       seed = 2, prior_only = TRUE)
head(edge_bayes_factors(fit, prior))
#>    type group1 group2 i j     prior posterior bayes_factor   log10_bf
#> 1 group      1     NA 1 2 0.2236667     0.515   3.68258205  0.5661524
#> 2 group      1     NA 1 3 0.2236667     0.415   2.46527463  0.3918653
#> 3 group      1     NA 2 3 0.2236667     0.800  13.74791098  1.1382367
#> 4 group      1     NA 1 4 0.2236667     0.080   0.31024463 -0.5082957
#> 5 group      1     NA 2 4 0.2236667     0.000   0.00865112 -2.0629277
#> 6 group      1     NA 3 4 0.2236667     0.035   0.13446120 -0.8714030
```
