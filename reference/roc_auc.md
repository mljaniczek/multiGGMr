# ROC curve and AUC from continuous scores

Computes the receiver operating characteristic (ROC) curve and area
under the curve (AUC) by comparing a continuous score matrix (e.g., PIP)
against a binary ground truth adjacency matrix. Only the upper triangle
is used.

## Usage

``` r
roc_auc(score_mat, truth_mat)
```

## Arguments

- score_mat:

  Numeric matrix `[p, p]`; continuous edge scores (e.g., posterior
  inclusion probabilities from
  [`pip_edges`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md)).

- truth_mat:

  Numeric matrix `[p, p]`; binary ground truth adjacency (1 = edge
  present, 0 = absent).

## Value

A list with three components:

- `FPR`:

  Numeric vector; false positive rates at each threshold (starts at 0,
  ends at 1).

- `TPR`:

  Numeric vector; true positive rates at each threshold (starts at 0,
  ends at 1).

- `auc`:

  Numeric scalar; area under the ROC curve. Values near 1 indicate
  excellent discrimination.

## See also

[`plot_roc()`](https://mljaniczek.github.io/multiGGMr/reference/plot_roc.md),
[`confusion_at_threshold()`](https://mljaniczek.github.io/multiGGMr/reference/confusion_at_threshold.md),
[`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
pip <- pip_edges(fit)
roc <- roc_auc(pip[,,1], sim$adj_list[[1]])
cat("AUC:", roc$auc, "\n")
#> AUC: 0.8564103 
```
