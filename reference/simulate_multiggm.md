# Simulate data from multiple Gaussian graphical models

Generates K precision matrices with controlled shared structure and
draws multivariate normal data from each. Follows the simulation design
in Peterson et al. (2015, JASA) Section 5.1. Group 1 uses the base graph
directly; groups 2 through K perturb the base graph by randomly flipping
edges with probability `perturb_prob`.

## Usage

``` r
simulate_multiggm(
  K = 2,
  p = 20,
  n = 100,
  graph_type = c("band", "random", "hub"),
  edge_prob = 0.1,
  perturb_prob = 0.1,
  signal = c(0.3, 0.6),
  seed = NULL
)
```

## Arguments

- K:

  Integer; number of sample groups. Default 2.

- p:

  Integer; number of variables (nodes). Default 20.

- n:

  Integer (scalar or length-K vector); sample size per group. If scalar,
  all groups have the same sample size. Default 100.

- graph_type:

  Character; type of base graph:

  - `band`: AR(2)-like banded structure: edges between all nodes within
    distance 2 (i.e., (i, i+1) and (i, i+2)). Produces approximately
    `2p - 3` edges. Matches Peterson et al. Section 5.1.

  - `random`: Erdos-Renyi random graph where each edge is included
    independently with probability `edge_prob`.

  - `hub`: Star/hub graph with `floor(p/5)` hub nodes, each connected to
    approximately 40\\

- edge_prob:

  Numeric; probability of each edge in the base graph (used only for
  `graph_type = "random"`). Default 0.1.

- perturb_prob:

  Numeric; for groups `k = 2, ..., K`, the probability that each edge is
  flipped (added if absent, removed if present) relative to the base
  graph. Controls how different the groups are. Set to 0 for identical
  graphs. Default 0.1.

- signal:

  Numeric vector of length 2; magnitude range `c(lo, hi)` for
  off-diagonal precision entries where edges exist. Signs are chosen
  randomly. Default `c(0.3, 0.6)`.

- seed:

  Optional integer random seed for reproducibility.

## Value

A list with components:

- `Omega_list`: List of K true precision matrices (each p x p, symmetric
  positive definite). Off-diagonal entries are non-zero only where edges
  exist, after row-normalization to ensure positive definiteness.

- `adj_list`: List of K true binary adjacency matrices (each p x p,
  0/1). `adj_list[[k]][i,j] = 1` if edge (i,j) exists in group k.

- `data_list`: List of K data matrices (each n_k x p), drawn from \\N(0,
  \Omega_k^{-1})\\ and column-centered.

- `S_list`: List of K cross-product matrices (each p x p), where
  `S_list[[k]] = t(X_k) %*% X_k` after centering.

- `n_vec`: Integer vector of length K with sample sizes.

- `K`: Integer; number of groups.

- `p`: Integer; number of variables.

## See also

[`multiggm_mcmc()`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 10, n = 100, seed = 42)
str(sim, max.level = 1)
#> List of 7
#>  $ Omega_list:List of 2
#>  $ adj_list  :List of 2
#>  $ data_list :List of 2
#>  $ S_list    :List of 2
#>  $ n_vec     : int [1:2] 100 100
#>  $ K         : num 2
#>  $ p         : num 10

# Check true edge counts
sapply(sim$adj_list, function(a) sum(a[upper.tri(a)]))
#> [1] 17 15
```
