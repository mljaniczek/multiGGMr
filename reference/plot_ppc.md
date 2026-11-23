# Plot posterior predictive checks

Plot posterior predictive checks

## Usage

``` r
plot_ppc(x, type = c("global", "edges"))
```

## Arguments

- x:

  A `"multiggm_ppc"` object from
  [`posterior_predictive_check()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_predictive_check.md).

- type:

  `"global"` (default): for each global statistic, the distribution over
  draws of \\T(y^{rep}) - T(y)\\; the ppp is the share above zero.
  `"edges"`: histograms of the edge-level ppp values per group and pair,
  which should not pile up at 0 or 1.

## Value

A `ggplot` object.

## See also

[`posterior_predictive_check()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_predictive_check.md)
