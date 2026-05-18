# Plot ROC curve

Creates a ROC curve plot from the output of
[`roc_auc`](https://mljaniczek.github.io/multiGGMr/reference/roc_auc.md),
with AUC annotation.

## Usage

``` r
plot_roc(roc_obj, main = "ROC Curve")
```

## Arguments

- roc_obj:

  A list as returned by
  [`roc_auc`](https://mljaniczek.github.io/multiGGMr/reference/roc_auc.md),
  with components `FPR`, `TPR`, and `auc`.

- main:

  Character; plot title. Default `"ROC Curve"`.

## Value

A `ggplot` object.

## See also

[`roc_auc()`](https://mljaniczek.github.io/multiGGMr/reference/roc_auc.md),
[`plot_recovery()`](https://mljaniczek.github.io/multiGGMr/reference/plot_recovery.md)

## Examples

``` r
sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
pip <- pip_edges(fit)
roc <- roc_auc(pip[,,1], sim$adj_list[[1]])
plot_roc(roc)

```
