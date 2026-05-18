# Simulate multi-platform data from multiple Gaussian graphical models

Generates S platforms of data, each with K groups sharing structure
across groups (controlled by `perturb_prob`) and across platforms
(controlled by `platform_perturb_prob`). Follows Shaddox et al. (2020)
simulation design. Platform 1 uses the base graph; platforms 2..S
perturb it.

## Usage

``` r
simulate_multiggm_platform(
  K = 2,
  S = 2,
  p_vec = c(20, 15),
  n = 100,
  graph_type = c("band", "random", "hub"),
  perturb_prob = 0.1,
  platform_perturb_prob = 0.2,
  signal = c(0.3, 0.6),
  seed = NULL
)
```

## Arguments

- K:

  Integer; number of sample groups. Default 2.

- S:

  Integer; number of platforms (data types). Default 2.

- p_vec:

  Integer vector of length S; number of variables per platform. Default
  `c(20, 15)`.

- n:

  Integer (scalar or length-K vector); sample size per group (same
  across platforms). Default 100.

- graph_type:

  Character; type of base graph. See
  [`simulate_multiggm`](https://mljaniczek.github.io/multiGGMr/reference/simulate_multiggm.md).

- perturb_prob:

  Numeric; probability of edge flip between groups within a platform.
  Default 0.1.

- platform_perturb_prob:

  Numeric; probability of edge flip between platforms. Default 0.2.

- signal:

  Numeric vector of length 2; magnitude range for off-diagonal precision
  entries. Default `c(0.3, 0.6)`.

- seed:

  Optional integer random seed.

## Value

A list with components:

- `platform_data`: List of S elements, each suitable as input to
  `multiggm_mcmc(method = "ssvs_platform", platform_data = ...)`. Each
  element is a list with `data_list`, `S_list`, `n_vec`, `Omega_list`,
  and `adj_list`.

- `K`, `S`, `p_vec`: Dimensions.

## See also

[`simulate_multiggm()`](https://mljaniczek.github.io/multiGGMr/reference/simulate_multiggm.md),
[`multiggm_mcmc()`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md)

## Examples

``` r
sim <- simulate_multiggm_platform(K = 2, S = 2, p_vec = c(10, 8),
                                   n = 80, seed = 42)
str(sim, max.level = 2)
#> List of 4
#>  $ platform_data:List of 2
#>   ..$ :List of 5
#>   ..$ :List of 5
#>  $ K            : num 2
#>  $ S            : num 2
#>  $ p_vec        : int [1:2] 10 8
```
