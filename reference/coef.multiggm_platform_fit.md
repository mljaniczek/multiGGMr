# Extract posterior mean precision matrices from a multi-platform fit

Extract posterior mean precision matrices from a multi-platform fit

## Usage

``` r
# S3 method for class 'multiggm_platform_fit'
coef(object, platform = NULL, ...)
```

## Arguments

- object:

  A `multiggm_platform_fit` object.

- platform:

  Integer; which platform to extract. If NULL (default), returns all
  platforms.

- ...:

  Ignored.

## Value

A list of K posterior mean precision matrices for the specified
platform, or a list of S such lists.
