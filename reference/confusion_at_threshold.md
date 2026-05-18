# Confusion metrics at a PIP threshold

Computes true/false positive/negative counts and rates by comparing a
continuous score matrix (e.g., PIP) against a binary ground truth
adjacency matrix at a given threshold. Only the upper triangle is used.

## Usage

``` r
confusion_at_threshold(score_mat, truth_mat, thr, upper_only = TRUE)
```

## Arguments

- score_mat:

  Numeric matrix `[p, p]`; continuous edge scores (e.g., posterior
  inclusion probabilities from
  [`pip_edges`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md)).

- truth_mat:

  Numeric matrix `[p, p]`; binary ground truth adjacency (1 = edge
  present, 0 = absent).

- thr:

  Numeric; threshold for edge selection. An edge is selected if
  `score_mat[i,j] >= thr`.

- upper_only:

  Logical; if `TRUE` (default), only the upper triangle is evaluated
  (appropriate for undirected graphs).

## Value

A named numeric vector with components:

- `TP`:

  True positives (correctly selected edges).

- `FP`:

  False positives (incorrectly selected non-edges).

- `TN`:

  True negatives (correctly excluded non-edges).

- `FN`:

  False negatives (missed true edges).

- `TPR`:

  True positive rate (sensitivity / recall).

- `FPR`:

  False positive rate (1 - specificity).

## See also

[`roc_auc()`](https://mljaniczek.github.io/multiGGMr/reference/roc_auc.md),
[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
pip <- pip_edges(fit)
confusion_at_threshold(pip[,,1], sim$adj_list[[1]], thr = 0.5)
#>         TP         FP         TN         FN        TPR        FPR 
#>  8.0000000  0.0000000 15.0000000  5.0000000  0.6153846  0.0000000 
```
