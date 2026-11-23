# Fit multiple Gaussian graphical models via MCMC

Implements Bayesian inference of multiple Gaussian graphical models
following Peterson, Stingo & Vannucci (2015, JASA). Simultaneously
estimates K group-specific precision matrices and their conditional
independence graphs, borrowing strength across groups through a Markov
random field (MRF) prior. Supports both raw data and pre-computed
cross-product matrices, with optional parallel multi-chain execution.

## Usage

``` r
multiggm_mcmc(
  data_list = NULL,
  S_list = NULL,
  n_vec = NULL,
  burnin = 5000,
  nsave = 1000,
  thin = 1,
  nchains = 1,
  parallel = FALSE,
  ncores = max(1L, parallel::detectCores() - 1L),
  seed = NULL,
  method = c("gwishart", "ssvs", "ssvs_platform", "mghs"),
  platform_data = NULL,
  hyper = NULL,
  engine = c("cpp", "R"),
  ...
)
```

## Arguments

- data_list:

  Optional list of K data matrices (each n_k x p). If provided, `S_list`
  and `n_vec` are computed automatically. Column-centers data
  internally. Exactly one of `data_list` or `S_list` must be given.

- S_list:

  Optional list of K cross-product matrices (each p x p), where
  `S_list[[k]] = t(X_k) %*% X_k` after column-centering.

- n_vec:

  Integer vector of length K with sample sizes. Required when `S_list`
  is provided; ignored when `data_list` is provided.

- burnin:

  Integer number of burn-in iterations (discarded). Default 5000.

- nsave:

  Integer number of posterior draws to save per chain (after burn-in and
  thinning). Default 1000.

- thin:

  Integer thinning interval; every `thin`-th post-burn-in iteration is
  saved. Default 1 (no thinning).

- nchains:

  Integer number of independent MCMC chains. Default 1.

- parallel:

  Logical; if `TRUE`, run chains in parallel using parallel. Default
  `FALSE`.

- ncores:

  Integer number of cores/workers when `parallel = TRUE`. Default is one
  less than the number of detected cores.

- seed:

  Optional integer random seed. If provided, chain `i` uses seed
  `seed + 1000 * (i - 1)` for reproducibility.

- method:

  Character; `"gwishart"` (default) uses the G-Wishart prior with the
  Wang-Li exchange algorithm (Peterson et al. 2015). `"ssvs"` uses a
  spike-and-slab normal prior on precision elements with column-wise
  Gibbs sampling (Shaddox et al. 2018), which scales to larger
  dimensions. `"ssvs_platform"` extends the SSVS approach to jointly
  model multiple data platforms (Shaddox et al. 2020, Biostatistics)
  with a hierarchical MRF coupling group similarity across platforms.
  Requires `platform_data`. `"mghs"` uses the multiple graphical
  horseshoe (Busatto & Stingo 2025): continuous horseshoe shrinkage on
  each precision matrix, with each edge's local scales correlated across
  groups. It has no edge-inclusion indicators; see Details.

- platform_data:

  For `method = "ssvs_platform"` only. A list of S elements (one per
  platform), each a list with components `data_list` or `S_list` and
  `n_vec`. Each platform may have different dimensions `p_s`, but must
  have the same number of groups K.

- hyper:

  Optional named list of hyperparameters. See Details.

- engine:

  Character; `"cpp"` (default) uses the high-performance C++ MCMC
  engine. `"R"` uses the pure R fallback (slower, mainly for debugging).
  Only applies when `method = "gwishart"`.

- ...:

  Additional arguments forwarded to the single-chain engine. For
  `method = "gwishart"`, `"ssvs"` and `"mghs"`: `save_draws` (default
  `TRUE`; `FALSE` keeps only the running summaries, see "Edge
  probabilities and memory") and `diff_deltas` (the thresholds `c` for
  the running probabilities \\P(\|x\| \> c)\\, \\P(\|D\| \> c)\\ and the
  union rule, default `c(0.01, 0.02, 0.03, 0.05, 0.10)`), `mcse_batches`
  (number of batches for the Monte Carlo SEs of edge probabilities,
  default 20; 0 turns them off), `waic` (`TRUE` keeps the
  per-observation sums for
  [`compute_waic()`](https://mljaniczek.github.io/multiGGMr/reference/compute_waic.md);
  needs `data_list`) and `prior_only` (`TRUE` ignores the data and
  samples the prior, for
  [`edge_bayes_factors()`](https://mljaniczek.github.io/multiGGMr/reference/edge_bayes_factors.md);
  not for mGHS). For `method = "mghs"` also `max_attempts` (restarts
  allowed when the chain gets stuck, default 3).

## Value

For a single chain (`nchains = 1`), an object of class `"multiggm_fit"`
with the following components:

- `K`:

  Integer; number of groups.

- `p`:

  Integer; number of variables (nodes).

- `C_save`:

  Numeric array `[p, p, K, nsave]`; posterior draws of the precision
  matrices \\\Omega_k\\. `C_save[,,k,s]` is the precision matrix for
  group `k` at saved iteration `s`.

- `Sig_save`:

  Numeric array `[p, p, K, nsave]`; posterior draws of the covariance
  matrices \\\Sigma_k = \Omega_k^{-1}\\. `Sig_save[,,k,s]` is the
  covariance matrix for group `k` at saved iteration `s`.

- `adj_save`:

  Integer array `[p, p, K, nsave]`; posterior draws of the adjacency
  matrices \\G_k\\. `adj_save[i,j,k,s]` is 1 if edge (i,j) is included
  in group `k` at iteration `s`, 0 otherwise.

- `Theta_save`:

  Numeric array `[K, K, nsave]`; posterior draws of the graph similarity
  parameters \\\theta\_{km}\\. `Theta_save[k,m,s]` measures similarity
  between graphs `k` and `m` at iteration `s`. Values \> 0 indicate the
  model borrows strength between those groups.

- `nu_save`:

  Numeric array `[p, p, nsave]`; posterior draws of the edge-specific
  log-odds parameters \\\nu\_{ij}\\. `nu_save[i,j,s]` controls the
  baseline edge inclusion probability for edge (i,j) across all groups
  at iteration `s`.

- `ar_gamma`:

  Numeric matrix `[K, K]`; acceptance rates for the between-model
  (spike-slab toggle) moves on \\\theta\_{km}\\. Upper triangle only.

- `ar_theta`:

  Numeric matrix `[K, K]`; acceptance rates for the within-model
  (slab-to-slab) moves on \\\theta\_{km}\\. Upper triangle only. Only
  meaningful when \\\theta\_{km} \> 0\\.

- `ar_nu`:

  Numeric matrix `[p, p]`; acceptance rates for the \\\nu\_{ij}\\ MH
  updates. Upper triangle only.

- `hyper`:

  Named list of the hyperparameters used.

- `call`:

  The matched function call.

For multiple chains (`nchains > 1`), an object of class
`"multiggm_fit_list"` with components:

- `chains`:

  List of `nchains` `multiggm_fit` objects.

- `call`:

  The matched function call.

- `K`:

  Integer; number of groups.

- `nchains`:

  Integer; number of chains.

## Details

**Method selection:** `method = "gwishart"` uses the G-Wishart prior on
precision matrices with the Wang-Li (2012) exchange algorithm. Accurate
but computationally expensive for \\p \> 30\\. `method = "ssvs"`
replaces this with a spike-and-slab normal prior, enabling column-wise
Gibbs sampling that scales to \\p \> 100\\. Both methods use the same
MRF prior for borrowing strength across groups.

The `hyper` list may contain any of the following:

**Common hyperparameters (both methods):**

- `a, b`:

  Beta(a, b) prior on \\q\_{ij} = \text{logit}^{-1}(\nu\_{ij})\\.
  Default a = 1, b = 4, giving prior edge probability \\\approx 0.20\\.

- `alpha, beta`:

  Gamma(shape = alpha, rate = beta) slab prior on \\\theta\_{km}\\.
  NOTE: `beta` is a RATE parameter. Default alpha = 2, beta = 5, giving
  mean = 0.4.

- `w`:

  Bernoulli prior probability that \\\gamma\_{km} = 1\\ (graph
  similarity active). Default 0.9 (strong prior belief that groups are
  related).

- `alpha_prop, beta_prop`:

  Gamma(shape = alpha_prop, rate = beta_prop) proposal for
  \\\theta\_{km}\\, a RATE like `beta` (all samplers). Default
  alpha_prop = 2, beta_prop = 5 for G-Wishart; alpha_prop = 1, beta_prop
  = 1 for SSVS.

**G-Wishart-specific hyperparameters (`method = "gwishart"`):**

- `b_prior`:

  G-Wishart degrees of freedom (default 3).

- `D_prior`:

  G-Wishart scale matrix, p x p (default `diag(p)`).

**SSVS-specific hyperparameters (`method = "ssvs"` or
`"ssvs_platform"`):**

- `v0`:

  Spike variance for off-diagonal precision elements (default
  `0.02^2 = 0.0004`).

- `v1`:

  Slab variance for off-diagonal precision elements (default
  `50^2 * v0 = 1.0`).

- `lambda`:

  Exponential rate hyperparameter for diagonal elements (default 1).

**mGHS hyperparameters (`method = "mghs"`):**

- `t_mean, t_var`:

  Mean and variance of the Beta prior on the cut-model threshold.
  Default 0.6 and 0.05, the values used by the authors' simulation
  scripts. The cut model does not feed back into the precision matrices.

- `R_eta`:

  Shape \\\eta\\ of an LKJ prior on the correlation of edge scales
  across groups, \\\pi(R) \propto \|R\|^{\eta - 1}\\. Default 1, the
  paper's uniform prior (the sampler is then unchanged). With few
  observations per group, the posterior of \\R\\ under the uniform prior
  piles up at \\\pm 1\\, because the per-group local scales can align
  the scaled edges of the groups and null edges carry no sign
  information; the chain then fails with "global-scale parameter above
  1e6" and an error message naming `R_eta`. In that case use \\\eta \>
  1\\: on a TCGA example (K = 2, n = 100, p = 28) \\\eta = 5\\ and
  \\\eta = 10\\ ran without failure where \\\eta \le 3\\ failed, and
  partial correlations and differential edge probabilities changed
  little between them, while the posterior of \\R\\ itself moved with
  \\\eta\\. Treat \\R\\ as a nuisance parameter in that regime. For K =
  2, \\R\_{12}\\ has a Beta(\\\eta\\, \\\eta\\) distribution stretched
  to (-1, 1).

**mGHS sampler (`method = "mghs"`):** the horseshoe prior is continuous,
so the sampler has no edge indicators: `adj_save` and `Sig_save` are
`NULL`, and
[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md),
[`pip_diff_edge()`](https://mljaniczek.github.io/multiGGMr/reference/pip_diff_edge.md),
[`edge_counts()`](https://mljaniczek.github.io/multiGGMr/reference/edge_counts.md),
`plot(type = "pip")` and
[`plot_recovery()`](https://mljaniczek.github.io/multiGGMr/reference/plot_recovery.md)
stop with an error. Select differential edges with
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)
(`rule = "magnitude"`). Do not call differences by whether a credible
interval excludes zero: under a continuous prior that rule is badly
anti-conservative (see
[`diff_summary()`](https://mljaniczek.github.io/multiGGMr/reference/diff_summary.md)).
The fit also contains `R_save` (`[K, K, nsave]`, correlation of edge
scales across groups), `tau_save` (`[K, nsave]`, global shrinkage),
`C_mean` and `pcor_mean` (posterior means over every post-burn-in
iteration), `diff_acc` (running sums of group differences), `cut_z_mean`
and `ta_save` (the authors' cut-model edge selection, kept for reference
only) and `attempts`. Storage is \\8 p^2 K\\ `nsave` bytes per chain
(160 MB at p = 100, K = 2, nsave = 1000); use `thin` to reduce `nsave`,
or `save_draws = FALSE`, in which case
[`diff_summary()`](https://mljaniczek.github.io/multiGGMr/reference/diff_summary.md),
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) use the running
summaries and draw-based functions stop with an error. A chain whose
shrinkage parameters become non-finite is restarted from scratch, up to
`max_attempts` times.

**Edge probabilities and memory.** Every sampler except
`"ssvs_platform"` keeps running sums over every post-burn-in iteration
(`edge_acc`) from which
[`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md),
[`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md)
and
[`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)
compute per-group and per-pair posterior edge probabilities, and
`C_mean` / `pcor_mean` (posterior means). They do not draw random
numbers, so a chain is the same with or without them. Their size grows
with the number of edges \\p(p-1)/2\\, the number of group pairs and the
number of thresholds in `diff_deltas`: with K = 2 and the default 5
thresholds, about 5 MB at p = 100 and 30 MB at p = 250, half of which
are the batch counts for Monte Carlo SEs (`mcse_batches = 0` drops
them). Saved draws take \\8 p^2 K\\ `nsave` bytes for the precision
matrices alone (160 MB at p = 100, K = 2, `nsave` = 1000); G-Wishart and
SSVS also save covariance, edge-indicator and \\\nu\\ draws of similar
size. With `save_draws = FALSE` only the running sums and the small
traces are kept; thresholds `c` must then come from `diff_deltas`, so
choose a fine grid if you want to pick `c` after fitting.

**Diagnostics.** Every sampler except `"ssvs_platform"` also keeps,
without drawing random numbers:

- `trace`: per post-burn-in iteration, the number of edges per group,
  the number of differential edges per pair of groups (structural,
  \\\|D\| \> c\\ and union, for each `c`) and the log-likelihood of each
  group, used by
  [`mcmc_diagnostics()`](https://mljaniczek.github.io/multiGGMr/reference/mcmc_diagnostics.md),
  [`plot_trace()`](https://mljaniczek.github.io/multiGGMr/reference/plot_trace.md)
  and
  [`diff_count_posterior()`](https://mljaniczek.github.io/multiGGMr/reference/diff_count_posterior.md);

- batch counts in `edge_acc$batch` for the Monte Carlo SE of every edge
  probability (`mcse_*` columns of
  [`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md)),
  which roughly double the size of `edge_acc`;

- node-level sums in `edge_acc$node` for
  [`node_diff_summary()`](https://mljaniczek.github.io/multiGGMr/reference/node_diff_summary.md);

- with `waic = TRUE`, `waic_acc` for
  [`compute_waic()`](https://mljaniczek.github.io/multiGGMr/reference/compute_waic.md):
  \\3 n_k\\ numbers per group, at a cost of \\O(n_k p^2)\\ per
  iteration.

The log-likelihood needs one Cholesky factorization per group and
iteration (\\O(p^3)\\), small next to a sampler sweep.

**Platform-level hyperparameters (`method = "ssvs_platform"` only):**

- `eta, kappa`:

  Gamma(shape = eta, scale = kappa) slab prior on platform similarity
  \\\Phi\_{st}\\. Default eta = 4, kappa = 5.

- `eta_prop, kappa_prop`:

  Gamma(shape = eta_prop, scale = kappa_prop) proposal for
  \\\Phi\_{st}\\, a SCALE like `kappa`. Default eta_prop = 1, kappa_prop
  = 1.

- `d, f`:

  Beta(d, f) prior on platform-level edge log-odds \\w\_{km}\\. Default
  d = 1, f = 19.

- `d_prop, f_prop`:

  Beta proposal for \\w\_{km}\\. Default d_prop = 2, f_prop = 4.

- `u_prior`:

  Bernoulli prior probability \\P(\zeta\_{st} = 1)\\. Default 0.1.

## References

Peterson, C.B., Stingo, F.C. & Vannucci, M. (2015). Bayesian inference
of multiple Gaussian graphical models. *Journal of the American
Statistical Association*, 110(509), 159-174.

Shaddox, E., Stingo, F.C., Peterson, C.B., et al. (2018). A Bayesian
approach for learning gene networks underlying disease severity in COPD.
*Statistics in Biosciences*, 10(1), 59-85.

Shaddox, E., Peterson, C.B., Stingo, F.C., et al. (2020). Bayesian
inference of networks across multiple sample groups and data types.
*Biostatistics*, 21(3), 561-576.

Busatto, C. & Stingo, F.C. (2025). Inference of multiple
high-dimensional networks with the Graphical Horseshoe prior. *Journal
of Computational and Graphical Statistics*, 1-12.
[doi:10.1080/10618600.2024.2435480](https://doi.org/10.1080/10618600.2024.2435480)

Busatto, C. & Stingo, F.C. (2025). The multiple graphical horseshoe: an
application to breast cancer data. In *SIS 2025*, 435-440. Springer.
[doi:10.1007/978-3-031-96736-8_72](https://doi.org/10.1007/978-3-031-96736-8_72)

## See also

[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md),
[`posterior_precision()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_precision.md),
[`posterior_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md),
[`summary.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/summary.multiggm_fit.md),
[`coef.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/coef.multiggm_fit.md),
[`fitted.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/fitted.multiggm_fit.md)

## Examples

``` r
# Short runs to keep the example fast; real analyses need far more
# iterations (see mcmc_diagnostics()).
sim <- simulate_multiggm(K = 2, p = 8, n = 100, seed = 42)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
summary(fit)
#> multiGGM MCMC Summary
#> =====================
#> Groups (K): 2   |  Nodes (p): 8   |  Posterior draws: 100 
#> 
#> Acceptance rates:
#>   gamma (edge toggle): 3.0% 
#>   theta (within-model): 44.3% 
#>   nu (edge log-odds): 56.4% 
#> 
#> Selected edges (PIP >= 0.5 ):
#>   Group 1 : 8 edges
#>   Group 2 : 9 edges
#> 
#> Graph similarity (theta):
#>   theta[1,2]: mean = 0.574, P(nonzero) = 99.0%
#> 

# Multiple chains (parallel on Unix/macOS, sequential on Windows)
fit_mc <- multiggm_mcmc(
  data_list = sim$data_list, burnin = 200, nsave = 100,
  nchains = 2, parallel = TRUE, ncores = 2, seed = 42
)
fit_mc                          # prints chain count and K
#> <multiggm_fit_list>
#>   Chains: 2 
#>   K groups: 2 
summary(fit_mc$chains[[1]])     # summarize first chain
#> multiGGM MCMC Summary
#> =====================
#> Groups (K): 2   |  Nodes (p): 8   |  Posterior draws: 100 
#> 
#> Acceptance rates:
#>   gamma (edge toggle): 2.3% 
#>   theta (within-model): 32.8% 
#>   nu (edge log-odds): 55.8% 
#> 
#> Selected edges (PIP >= 0.5 ):
#>   Group 1 : 8 edges
#>   Group 2 : 11 edges
#> 
#> Graph similarity (theta):
#>   theta[1,2]: mean = 0.774, P(nonzero) = 100.0%
#> 
pip_edges(fit_mc, chain = 1)    # PIPs from chain 1
#> , , 1
#> 
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,] 0.00 1.00 0.99 0.21 0.17 0.10 0.11 0.70
#> [2,] 1.00 0.00 0.58 0.95 0.31 0.00 0.22 0.08
#> [3,] 0.99 0.58 0.00 0.33 0.33 0.11 0.33 0.00
#> [4,] 0.21 0.95 0.33 0.00 0.38 0.88 0.35 0.07
#> [5,] 0.17 0.31 0.33 0.38 0.00 0.15 0.63 0.06
#> [6,] 0.10 0.00 0.11 0.88 0.15 0.00 0.42 0.77
#> [7,] 0.11 0.22 0.33 0.35 0.63 0.42 0.00 0.48
#> [8,] 0.70 0.08 0.00 0.07 0.06 0.77 0.48 0.00
#> 
#> , , 2
#> 
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,] 0.00 0.42 1.00 0.07 0.13 0.04 0.12 0.79
#> [2,] 0.42 0.00 0.68 1.00 0.23 0.00 0.08 0.10
#> [3,] 1.00 0.68 0.00 0.59 0.99 0.09 0.53 0.00
#> [4,] 0.07 1.00 0.59 0.00 0.28 1.00 0.20 0.07
#> [5,] 0.13 0.23 0.99 0.28 0.00 0.12 0.34 0.27
#> [6,] 0.04 0.00 0.09 1.00 0.12 0.00 0.98 0.96
#> [7,] 0.12 0.08 0.53 0.20 0.34 0.98 0.00 0.57
#> [8,] 0.79 0.10 0.00 0.07 0.27 0.96 0.57 0.00
#> 
pip_edges(fit_mc, chain = 2)    # PIPs from chain 2
#> , , 1
#> 
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,] 0.00 1.00 0.96 0.14 0.07 0.13 0.02 0.86
#> [2,] 1.00 0.00 0.49 0.63 0.53 0.07 0.25 0.16
#> [3,] 0.96 0.49 0.00 0.37 0.23 0.08 0.31 0.12
#> [4,] 0.14 0.63 0.37 0.00 0.50 0.85 0.49 0.11
#> [5,] 0.07 0.53 0.23 0.50 0.00 0.06 0.26 0.09
#> [6,] 0.13 0.07 0.08 0.85 0.06 0.00 0.33 0.61
#> [7,] 0.02 0.25 0.31 0.49 0.26 0.33 0.00 0.38
#> [8,] 0.86 0.16 0.12 0.11 0.09 0.61 0.38 0.00
#> 
#> , , 2
#> 
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,] 0.00 0.30 1.00 0.15 0.02 0.12 0.07 0.84
#> [2,] 0.30 0.00 0.36 1.00 0.51 0.08 0.14 0.09
#> [3,] 1.00 0.36 0.00 0.79 0.98 0.03 0.69 0.05
#> [4,] 0.15 1.00 0.79 0.00 0.31 1.00 0.33 0.12
#> [5,] 0.02 0.51 0.98 0.31 0.00 0.13 0.08 0.27
#> [6,] 0.12 0.08 0.03 1.00 0.13 0.00 0.96 0.91
#> [7,] 0.07 0.14 0.69 0.33 0.08 0.96 0.00 0.44
#> [8,] 0.84 0.09 0.05 0.12 0.27 0.91 0.44 0.00
#> 
```
