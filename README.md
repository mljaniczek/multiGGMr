
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multiGGM

<!-- badges: start -->

<!-- badges: end -->

The **multiGGM** package implements the Bayesian method of Peterson,
Stingo & Vannucci (2015, *JASA*) for joint inference of multiple
Gaussian graphical models (GGMs). Given data from $K$ groups, the method
simultaneously estimates group-specific precision matrices $\Omega_k$
and their associated conditional independence graphs $G_k$, while
borrowing strength across groups through a Markov random field (MRF)
prior. It also implements two follow-up versions of the original
algorithm, a faster version using continuous spike-and-slab priors
called scale-multiGGM (Shaddox, Peterson, Vanucci, 2018, *Statistics in
Biosciences*), and multi-platform multiGGM (Shaddox et al. 2020) which
allows you to incorporate multi-modal data and jointly estimate
networks.

## Installation

You can install the development version of multiGGM from
[GitHub](https://github.com/mljaniczek/multiGGMr) with:

``` r
# install.packages("pak")
pak::pak("mljaniczek/multiGGMr")
```

## Example

This is a basic example simulating data from 2 groups and then running
the scale-multiGGM. For more detailed examples see vignettes.

``` r
library(multiGGM)

# first simulate some data
sim <- simulate_multiggm(
  K = 2, p = 20, n = 100,
  graph_type = "hub",
  perturb_prob = 0.1,
  seed = 42
)

# then run multiGGM, specifying which method

fit <- multiggm_mcmc(
  data_list = sim$data_list,
  method    = "ssvs", # this is for the scalable version
  burnin    = 200,
  nsave     = 100,
  seed      = 123
)
```

## References

Peterson, C.B., Stingo, F.C. & Vannucci, M. (2015). Bayesian inference
of multiple Gaussian graphical models. *Journal of the American
Statistical Association*, 110(509), 159–174.

Shaddox, E., Stingo, F.C., Peterson, C.B., et al. (2018). A Bayesian
approach for learning gene networks underlying disease severity in COPD.
*Statistics in Biosciences*, 10(1), 59–85.

Shaddox, E., Peterson, C.B., Stingo, F.C., et al. (2020). Bayesian
inference of networks across multiple sample groups and data types.
*Biostatistics*, 21 (3) 561–-576.

Wang, H. (2015). Scaling it up: stochastic search structure learning in
graphical models. *Bayesian Analysis*, 10(2), 351–377.

Wang, H. & Li, S.Z. (2012). Efficient Gaussian graphical model
determination under G-Wishart prior distributions. *Electronic Journal
of Statistics*, 6, 168–198.
