# Calibration of Bayesian FDR selection against a known truth

For simulation studies: selects edges by Bayesian FDR at each target
level `alpha` and compares the posterior expected FDR with the achieved
false discovery proportion and power under the true precision matrices.
Well-calibrated edge probabilities give achieved FDR close to the
expected FDR.

## Usage

``` r
fdr_calibration(
  fit,
  truth,
  groups = c(1, 2),
  group = NULL,
  rule = c("magnitude", "structural", "union", "cut"),
  c = NULL,
  truth_c = c,
  scale = c("pcor", "precision"),
  alpha = seq(0.01, 0.2, by = 0.01),
  tol = 1e-08,
  chain = NULL
)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list` object.

- truth:

  A list of the K true precision matrices, or the output of
  [`simulate_multiggm()`](https://mljaniczek.github.io/multiGGMr/reference/simulate_multiggm.md).

- groups:

  Integer vector of length 2, the groups `c(k, m)` to compare. Default
  `c(1, 2)`.

- group:

  A single group to calibrate within-group selection
  ([`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md));
  `NULL` (default) calibrates differential selection between `groups`.

- rule:

  Selection rule as in
  [`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)
  or
  [`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md).

- c:

  Threshold for the magnitude and union rules.

- truth_c:

  Threshold defining a true difference (or true edge) for the magnitude
  rule, on `scale`. Default `c`: a true difference is \\\|D\_{true}\| \>
  c\\. Structural truth is a difference in the nonzero pattern (or a
  nonzero entry); union truth is either.

- scale:

  `"pcor"` (default) for partial correlations or `"precision"` for
  precision-matrix entries.

- alpha:

  Vector of target FDR levels.

- tol:

  Entries with absolute value at most `tol` count as zero in the true
  precision matrices.

- chain:

  `NULL` (default) pools all chains; an integer vector selects chains of
  a `multiggm_fit_list`.

## Value

A data frame with `alpha`, `n_selected`, `expected_fdr` (posterior
expected FDR of the selection), `achieved_fdr` (false discovery
proportion), `power` and `n_true`. Attributes `rule`, `c`, `truth_c`.

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 100, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs",
                     burnin = 200, nsave = 200, save_draws = FALSE)
fdr_calibration(fit, sim, rule = "union", c = 0.1)
#>    alpha n_selected expected_fdr achieved_fdr     power n_true
#> 1   0.01          3  0.003333333    0.0000000 0.3333333      9
#> 2   0.02          4  0.015000000    0.2500000 0.3333333      9
#> 3   0.03          4  0.015000000    0.2500000 0.3333333      9
#> 4   0.04          5  0.034000000    0.4000000 0.3333333      9
#> 5   0.05          5  0.034000000    0.4000000 0.3333333      9
#> 6   0.06          5  0.034000000    0.4000000 0.3333333      9
#> 7   0.07          5  0.034000000    0.4000000 0.3333333      9
#> 8   0.08          5  0.034000000    0.4000000 0.3333333      9
#> 9   0.09          5  0.034000000    0.4000000 0.3333333      9
#> 10  0.10          5  0.034000000    0.4000000 0.3333333      9
#> 11  0.11          6  0.109166667    0.3333333 0.4444444      9
#> 12  0.12          6  0.109166667    0.3333333 0.4444444      9
#> 13  0.13          6  0.109166667    0.3333333 0.4444444      9
#> 14  0.14          6  0.109166667    0.3333333 0.4444444      9
#> 15  0.15          6  0.109166667    0.3333333 0.4444444      9
#> 16  0.16          6  0.109166667    0.3333333 0.4444444      9
#> 17  0.17          6  0.109166667    0.3333333 0.4444444      9
#> 18  0.18          7  0.172857143    0.2857143 0.5555556      9
#> 19  0.19          7  0.172857143    0.2857143 0.5555556      9
#> 20  0.20          7  0.172857143    0.2857143 0.5555556      9
```
