# Simulate multi-platform data from multiple Gaussian graphical models

Generates S platforms (data types) of data, each with the same K groups,
with graphs that share structure across platforms and across groups.

## Usage

``` r
simulate_multiggm_platform(
  K = 2,
  S = 2,
  p_vec = c(20, 15),
  n = 100,
  graph_type = c("band", "random", "hub"),
  edge_prob = 0.1,
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

- edge_prob:

  Numeric; edge probability of the base graph for
  `graph_type = "random"`. Default 0.1.

- perturb_prob:

  Numeric; probability of edge flip between groups within a platform.
  Default 0.1.

- platform_perturb_prob:

  Numeric; probability of edge flip between platform 1 and each other
  platform. Default 0.2.

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

- `adj_base`: The base graph on `max(p_vec)` nodes.

- `platform_adj`: List of S platform graphs (before group perturbation).

- `K`, `S`, `p_vec`: Dimensions.

## Details

The data-generating process is hierarchical:

1.  One base graph is drawn on `max(p_vec)` nodes (see `graph_type`).

2.  Platform `s` uses the leading `p_vec[s] x p_vec[s]` block of the
    base graph. Platform 1 uses it unchanged; platforms 2, ..., S flip
    each pair independently with probability `platform_perturb_prob`.
    Variable `j` therefore plays the same role in every platform that
    has it.

3.  Within each platform, group 1 uses the platform graph and groups 2,
    ..., K flip each pair with probability `perturb_prob`, exactly as in
    [`simulate_multiggm()`](https://mljaniczek.github.io/multiGGMr/reference/simulate_multiggm.md).
    Group perturbations are independent across platforms.

4.  Precision matrices and data are generated from each graph as in
    [`simulate_multiggm()`](https://mljaniczek.github.io/multiGGMr/reference/simulate_multiggm.md).

## See also

[`simulate_multiggm()`](https://mljaniczek.github.io/multiGGMr/reference/simulate_multiggm.md),
[`multiggm_mcmc()`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md)

## Examples

``` r
sim <- simulate_multiggm_platform(K = 2, S = 2, p_vec = c(10, 8),
                                   n = 80, seed = 42)
str(sim, max.level = 2)
#> List of 6
#>  $ platform_data:List of 2
#>   ..$ :List of 5
#>   ..$ :List of 5
#>  $ adj_base     : int [1:10, 1:10] 0 1 1 0 0 0 0 0 0 0 ...
#>  $ platform_adj :List of 2
#>   ..$ : int [1:10, 1:10] 0 1 1 0 0 0 0 0 0 0 ...
#>   ..$ : int [1:8, 1:8] 0 0 1 0 0 0 0 0 0 0 ...
#>  $ K            : num 2
#>  $ S            : num 2
#>  $ p_vec        : int [1:2] 10 8
```
