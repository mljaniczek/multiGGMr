# Internal multi-platform SSVS engine

Implements Shaddox et al. (2020) multi-platform SSVS with hierarchical
MRF coupling across S data platforms.

## Usage

``` r
.multiggm_mcmc_platform(
  platform_data,
  burnin,
  nsave,
  thin,
  nchains = 1L,
  parallel = FALSE,
  ncores = 1L,
  seed = NULL,
  hyper = NULL,
  disp = FALSE,
  ...
)
```
