# MCMC convergence diagnostics

Collects convergence diagnostics for scalar summaries of the chain, the
Monte Carlo error of edge probabilities, agreement between chains and
the samplers' acceptance rates in one table. Works with
`save_draws = FALSE`: the traces and batch counts it uses are always
kept.

## Usage

``` r
mcmc_diagnostics(
  fit,
  c = NULL,
  scale = c("pcor", "precision"),
  rhat_max = 1.01,
  ess_min = 400,
  mcse_max = 0.05
)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list`. Several chains are strongly
  recommended: R-hat from a single chain only compares its two halves.

- c:

  Threshold for the trace of the number of differential edges \\\|D\| \>
  c\\ and for the Monte Carlo error of \\P(\|D\| \> c)\\. Default: the
  largest value of `diff_deltas`.

- scale:

  Scale of \\D\\: `"pcor"` (default) or `"precision"`.

- rhat_max, ess_min, mcse_max:

  Thresholds for the flags: R-hat above `rhat_max` (default 1.01), bulk
  or tail ESS below `ess_min` (default 400) and edge-probability Monte
  Carlo SE above `mcse_max` (default 0.05).

## Value

An object of class `"multiggm_diagnostics"`, a list with

- `scalars`:

  Data frame with one row per scalar summary: `quantity`, `mean`, `sd`,
  `rhat`, `ess_bulk`, `ess_tail`, `mcse_mean` and `flag`.

- `edges`:

  Data frame with one row per edge probability: `quantity`, `max_mcse`,
  `median_mcse`, `n_mcse_above`, `max_chain_diff`.

- `acceptance`:

  Data frame of mean acceptance rates for the Metropolis-Hastings steps
  (G-Wishart and SSVS), with a flag outside (0.05, 0.95).

- `n_chains`, `n_iter`, `c`, `scale`, `thresholds`:

  Settings.

## Details

**Scalar traces.** For each chain the samplers record, at every
post-burn-in iteration, the number of edges per group, the number of
differential edges per pair of groups (structural, \\\|D\| \> c\\ and
union for every `c` in `diff_deltas`) and the Gaussian log-likelihood of
each group. These, together with \\\theta\\ (G-Wishart, SSVS) or the
edge-scale correlations \\R\\ and global scales \\\tau\\ (mGHS), are
summarized by rank-normalized split R-hat and bulk and tail effective
sample sizes (Vehtari et al. 2021), the same algorithm as the posterior
package. A whole-graph summary such as the number of edges is a
sensitive check: chains stuck in different parts of graph space disagree
on it.

**Edge probabilities.** R-hat is not informative for single 0/1 edge
indicators. Instead, `edges` reports the batch-means Monte Carlo SE of
every edge probability used for selection (largest and median over
edges, and how many exceed `mcse_max`), and, for several chains, the
largest difference in an edge probability between any two chains. An
edge whose probability is within a few Monte Carlo SEs of a selection
cutoff is not reliably selected or rejected; run the chain longer.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. & Buerkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved R-hat
for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718.

## See also

[`plot_trace()`](https://mljaniczek.github.io/multiGGMr/reference/plot_trace.md),
[`plot_chain_agreement()`](https://mljaniczek.github.io/multiGGMr/reference/plot_chain_agreement.md),
[`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs", burnin = 200,
                     nsave = 200, nchains = 2, seed = 1, save_draws = FALSE)
mcmc_diagnostics(fit)
#> <multiggm_diagnostics> ssvs | 2 chain(s) x 200 post-burn-in iterations
#> 
#> Scalar summaries (flag: R-hat > 1.01 or ESS < 400 ):
#>                quantity     mean    sd  rhat ess_bulk ess_tail mcse_mean flag
#>              n_edges[1]    5.410 1.051 1.154    9.199    8.678     0.348 rhat
#>               loglik[1] -503.899 3.504 1.080   18.663  104.940     0.899 rhat
#>              n_edges[2]    4.640 1.046 1.110   18.961   19.200     0.225 rhat
#>               loglik[2] -532.755 3.274 1.144   11.230   95.167     1.004 rhat
#>      n_diff_struct[1,2]    4.040 1.442 1.062   38.063   42.465     0.230 rhat
#>  n_diff_exceed_0.1[1,2]    6.068 1.232 1.051   24.034   54.545     0.256 rhat
#>   n_diff_union_0.1[1,2]    6.370 1.394 1.038   29.068   47.886     0.262 rhat
#>              theta[1,2]    0.448 0.293 1.054   28.597  165.368     0.060 rhat
#> 
#> Edge probabilities (Monte Carlo SE; n_mcse_above: edges with SE > 0.05 ):
#>             quantity max_mcse median_mcse n_mcse_above max_chain_diff
#>      group: p_struct    0.065       0.019            7          0.470
#>  group: p_exceed_0.1    0.061       0.012            5          0.460
#>       pair: p_struct    0.060       0.040            5          0.295
#>   pair: p_exceed_0.1    0.065       0.015            1          0.315
#>    pair: p_union_0.1    0.066       0.024            2          0.330
#> 
#> Acceptance rates:
#>                         step  rate flag
#>  gamma (between-model theta) 0.090     
#>         theta (within-model) 0.462     
#>                           nu 0.542     
#> 
#> 13 flag(s); see above.
```
