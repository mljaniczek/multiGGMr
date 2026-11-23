# One posterior edge probability as a p x p matrix

Extracts a single quantity from
[`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md)
for one group or one pair of groups, as a symmetric matrix.

## Usage

``` r
edge_prob_matrix(
  fit,
  quantity = c("p_struct", "p_cut", "mean_weight", "mean", "sd", "p_sign", "p_exceed",
    "p_union", "p_up", "p_down", "mcse_struct", "mcse_cut", "mcse_exceed", "mcse_union",
    "mean_shrinkage", "p_shrinkage", "mcse_shrinkage"),
  group = NULL,
  pair = NULL,
  c = NULL,
  scale = c("precision", "pcor"),
  chain = NULL
)
```

## Arguments

- fit:

  A `multiggm_fit` or `multiggm_fit_list`.

- quantity:

  One of the columns of
  [`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md)
  without the threshold suffix: `"p_struct"`, `"p_cut"`,
  `"mean_weight"`, `"mean"`, `"sd"`, `"p_sign"`, `"p_exceed"`,
  `"p_union"`, `"p_up"`, `"p_down"`, `"mcse_struct"`, `"mcse_cut"`,
  `"mcse_exceed"` or `"mcse_union"`.

- group:

  A single group, for per-group quantities.

- pair:

  Two groups `c(s, t)`, for differential quantities. For `"mean"`,
  `"p_up"` and `"p_down"` the direction follows the order (\\D = x_s -
  x_t\\).

- c:

  A single threshold, required for the quantities with a threshold
  (`"p_exceed"`, `"p_union"`, `"p_up"`, `"p_down"`, `"mcse_exceed"`,
  `"mcse_union"`).

- scale:

  `"precision"` (default) or `"pcor"`: the scale of `mean`, `sd`,
  `p_exceed_*`, `p_union_*` and `p_sign`.

- chain:

  `NULL` (default) pools all chains; an integer selects chains.

## Value

A symmetric `p x p` matrix with zero diagonal.

## See also

[`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md),
[`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 100, nsave = 100)
edge_prob_matrix(fit, "p_union", pair = c(1, 2), c = 0.05)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,] 0.00 0.50 0.79 0.15 0.32 0.13
#> [2,] 0.50 0.00 1.00 0.75 0.27 0.90
#> [3,] 0.79 1.00 0.00 0.31 0.98 0.40
#> [4,] 0.15 0.75 0.31 0.00 1.00 0.90
#> [5,] 0.32 0.27 0.98 1.00 0.00 1.00
#> [6,] 0.13 0.90 0.40 0.90 1.00 0.00
```
