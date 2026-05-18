# Compare estimated graph to ground truth

Shows side-by-side heatmaps of the true adjacency matrix, PIP, and
thresholded estimate. Returns confusion metrics invisibly.

## Usage

``` r
plot_recovery(fit, true_adj, pip_threshold = 0.5, groups = NULL)
```

## Arguments

- fit:

  A `multiggm_fit` object returned by
  [`multiggm_mcmc`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md).

- true_adj:

  A list of K true adjacency matrices (0/1), or a single matrix
  (recycled for all groups).

- pip_threshold:

  Numeric; threshold for edge selection. Default 0.5.

- groups:

  Integer vector; which groups to plot. Default: all groups.

## Value

A named list of confusion metric vectors per group (invisible). Each
element is a named numeric vector as returned by
[`confusion_at_threshold`](https://mljaniczek.github.io/multiGGMr/reference/confusion_at_threshold.md),
with components `TP`, `FP`, `TN`, `FN`, `TPR`, `FPR`.

## See also

[`confusion_at_threshold()`](https://mljaniczek.github.io/multiGGMr/reference/confusion_at_threshold.md),
[`roc_auc()`](https://mljaniczek.github.io/multiGGMr/reference/roc_auc.md),
[`plot_roc()`](https://mljaniczek.github.io/multiGGMr/reference/plot_roc.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
cm <- plot_recovery(fit, sim$adj_list)

cm$Group_1  # confusion metrics for group 1
#>         TP         FP         TN         FN        TPR        FPR 
#>  8.0000000  0.0000000 15.0000000  5.0000000  0.6153846  0.0000000 
```
