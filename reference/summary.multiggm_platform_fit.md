# Summarize a multiggm_platform_fit object

Displays per-platform MCMC diagnostics and platform-level similarity
estimates.

## Usage

``` r
# S3 method for class 'multiggm_platform_fit'
summary(object, pip_threshold = 0.5, ...)
```

## Arguments

- object:

  A `multiggm_platform_fit` object.

- pip_threshold:

  Numeric threshold for counting selected edges. Default 0.5.

- ...:

  Ignored.

## Value

An object of class `"summary.multiggm_platform_fit"` (printed
invisibly).
