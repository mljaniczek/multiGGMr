# Extract posterior mean partial correlation matrices from a multi-platform fit

Extract posterior mean partial correlation matrices from a
multi-platform fit

## Usage

``` r
# S3 method for class 'multiggm_platform_fit'
fitted(object, platform = NULL, ...)
```

## Arguments

- object:

  A `multiggm_platform_fit` object.

- platform:

  Integer; which platform. If NULL, returns all.

- ...:

  Ignored.

## Value

A list of K posterior mean partial correlation matrices, or a list of S
such lists.
