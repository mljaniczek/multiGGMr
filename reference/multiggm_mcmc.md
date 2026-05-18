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
  method = c("gwishart", "ssvs", "ssvs_platform"),
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
  Requires `platform_data`.

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

  Additional arguments forwarded to the single-chain engine.

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
  \\\theta\_{km}\\. Default alpha_prop = 2, beta_prop = 5 for G-Wishart;
  alpha_prop = 1, beta_prop = 1 for SSVS.

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

**Platform-level hyperparameters (`method = "ssvs_platform"` only):**

- `eta, kappa`:

  Gamma(shape = eta, scale = kappa) slab prior on platform similarity
  \\\Phi\_{st}\\. Default eta = 4, kappa = 5.

- `eta_prop, kappa_prop`:

  Gamma proposal for \\\Phi\_{st}\\. Default eta_prop = 1, kappa_prop =
  1.

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

## See also

[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md),
[`posterior_precision()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_precision.md),
[`posterior_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md),
[`summary.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/summary.multiggm_fit.md),
[`coef.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/coef.multiggm_fit.md),
[`fitted.multiggm_fit()`](https://mljaniczek.github.io/multiGGMr/reference/fitted.multiggm_fit.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 10, n = 100, seed = 42)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 500, nsave = 200)
summary(fit)
#> multiGGM MCMC Summary
#> =====================
#> Groups (K): 2   |  Nodes (p): 10   |  Posterior draws: 200 
#> 
#> Acceptance rates:
#>   gamma (edge toggle): 7.9% 
#>   theta (within-model): 57.8% 
#>   nu (edge log-odds): 51.0% 
#> 
#> Selected edges (PIP >= 0.5 ):
#>   Group 1 : 10 edges
#>   Group 2 : 9 edges
#> 
#> Graph similarity (theta):
#>   theta[1,2]: mean = 0.416, P(nonzero) = 94.0%
#> 

# Multiple chains (parallel on Unix/macOS, sequential on Windows)
fit_mc <- multiggm_mcmc(
  data_list = sim$data_list, burnin = 500, nsave = 200,
  nchains = 2, parallel = TRUE, seed = 42
)
fit_mc                          # prints chain count and K
#> <multiggm_fit_list>
#>   Chains: 2 
#>   K groups: 2 
summary(fit_mc$chains[[1]])     # summarize first chain
#> multiGGM MCMC Summary
#> =====================
#> Groups (K): 2   |  Nodes (p): 10   |  Posterior draws: 200 
#> 
#> Acceptance rates:
#>   gamma (edge toggle): 5.0% 
#>   theta (within-model): 50.3% 
#>   nu (edge log-odds): 49.9% 
#> 
#> Selected edges (PIP >= 0.5 ):
#>   Group 1 : 10 edges
#>   Group 2 : 10 edges
#> 
#> Graph similarity (theta):
#>   theta[1,2]: mean = 0.603, P(nonzero) = 98.5%
#> 
pip_edges(fit_mc, chain = 1)    # PIPs from chain 1
#> , , 1
#> 
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
#>  [1,] 0.000 0.840 0.065 0.355 0.045 0.050 0.065 0.060 0.150 0.090
#>  [2,] 0.840 0.000 0.715 0.855 0.085 0.085 0.100 0.110 0.035 0.080
#>  [3,] 0.065 0.715 0.000 0.230 0.295 0.090 0.365 0.140 0.060 0.025
#>  [4,] 0.355 0.855 0.230 0.000 0.840 1.000 0.150 0.005 0.080 0.090
#>  [5,] 0.045 0.085 0.295 0.840 0.000 0.620 0.230 0.050 0.115 0.170
#>  [6,] 0.050 0.085 0.090 1.000 0.620 0.000 0.200 0.635 0.105 0.050
#>  [7,] 0.065 0.100 0.365 0.150 0.230 0.200 0.000 0.215 0.475 0.130
#>  [8,] 0.060 0.110 0.140 0.005 0.050 0.635 0.215 0.000 0.625 0.980
#>  [9,] 0.150 0.035 0.060 0.080 0.115 0.105 0.475 0.625 0.000 1.000
#> [10,] 0.090 0.080 0.025 0.090 0.170 0.050 0.130 0.980 1.000 0.000
#> 
#> , , 2
#> 
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
#>  [1,] 0.000 0.655 0.075 0.200 0.050 0.045 0.055 0.050 0.130 0.255
#>  [2,] 0.655 0.000 0.615 1.000 0.160 0.035 0.055 0.090 0.015 0.105
#>  [3,] 0.075 0.615 0.000 0.085 0.190 0.120 0.110 0.200 0.055 0.010
#>  [4,] 0.200 1.000 0.085 0.000 0.355 0.715 0.085 0.005 0.095 0.080
#>  [5,] 0.050 0.160 0.190 0.355 0.000 1.000 0.375 0.125 0.105 0.555
#>  [6,] 0.045 0.035 0.120 0.715 1.000 0.000 0.245 0.285 0.090 0.145
#>  [7,] 0.055 0.055 0.110 0.085 0.375 0.245 0.000 0.655 0.975 0.205
#>  [8,] 0.050 0.090 0.200 0.005 0.125 0.285 0.655 0.000 0.255 0.965
#>  [9,] 0.130 0.015 0.055 0.095 0.105 0.090 0.975 0.255 0.000 1.000
#> [10,] 0.255 0.105 0.010 0.080 0.555 0.145 0.205 0.965 1.000 0.000
#> 
pip_edges(fit_mc, chain = 2)    # PIPs from chain 2
#> , , 1
#> 
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
#>  [1,] 0.000 0.880 0.300 0.255 0.055 0.020 0.185 0.085 0.050 0.095
#>  [2,] 0.880 0.000 0.795 0.815 0.025 0.085 0.025 0.040 0.085 0.055
#>  [3,] 0.300 0.795 0.000 0.215 0.155 0.095 0.090 0.160 0.100 0.080
#>  [4,] 0.255 0.815 0.215 0.000 0.895 1.000 0.105 0.065 0.030 0.055
#>  [5,] 0.055 0.025 0.155 0.895 0.000 0.790 0.145 0.120 0.045 0.075
#>  [6,] 0.020 0.085 0.095 1.000 0.790 0.000 0.280 0.600 0.035 0.150
#>  [7,] 0.185 0.025 0.090 0.105 0.145 0.280 0.000 0.200 0.400 0.105
#>  [8,] 0.085 0.040 0.160 0.065 0.120 0.600 0.200 0.000 0.525 0.960
#>  [9,] 0.050 0.085 0.100 0.030 0.045 0.035 0.400 0.525 0.000 0.980
#> [10,] 0.095 0.055 0.080 0.055 0.075 0.150 0.105 0.960 0.980 0.000
#> 
#> , , 2
#> 
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
#>  [1,] 0.000 0.770 0.205 0.185 0.070 0.095 0.060 0.080 0.090 0.170
#>  [2,] 0.770 0.000 0.700 1.000 0.100 0.045 0.050 0.075 0.100 0.080
#>  [3,] 0.205 0.700 0.000 0.080 0.110 0.105 0.025 0.205 0.105 0.100
#>  [4,] 0.185 1.000 0.080 0.000 0.340 0.760 0.060 0.015 0.040 0.040
#>  [5,] 0.070 0.100 0.110 0.340 0.000 1.000 0.295 0.195 0.030 0.350
#>  [6,] 0.095 0.045 0.105 0.760 1.000 0.000 0.290 0.250 0.045 0.295
#>  [7,] 0.060 0.050 0.025 0.060 0.295 0.290 0.000 0.535 0.985 0.130
#>  [8,] 0.080 0.075 0.205 0.015 0.195 0.250 0.535 0.000 0.165 0.990
#>  [9,] 0.090 0.100 0.105 0.040 0.030 0.045 0.985 0.165 0.000 1.000
#> [10,] 0.170 0.080 0.100 0.040 0.350 0.295 0.130 0.990 1.000 0.000
#> 
```
