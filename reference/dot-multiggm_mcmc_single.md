# Internal single-chain MCMC engine

Thin R wrapper that calls the C++ MCMC engine. Falls back to pure R if
`engine = "R"`.

## Usage

``` r
.multiggm_mcmc_single(
  S_list,
  n_vec,
  burnin,
  nsave,
  thin,
  chain_id = 1L,
  Theta_init = NULL,
  nu_init = NULL,
  C_init = NULL,
  disp = FALSE,
  method = "gwishart",
  hyper = NULL,
  engine = "cpp",
  save_draws = TRUE,
  diff_deltas = .default_deltas,
  mcse_batches = 20L,
  waic = FALSE,
  prior_only = FALSE,
  X_list = NULL,
  ...
)
```
