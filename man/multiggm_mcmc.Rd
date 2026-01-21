\name{multiggm_mcmc}
\alias{multiggm_mcmc}
\title{Fit multi-graph Gaussian graphical model via MCMC}
\usage{
multiggm_mcmc(S_list, n_vec, burnin = 5000, nsave = 1000, thin = 1, nchains = 1,
  parallel = FALSE, ncores = max(1L, parallel::detectCores() - 1L), seed = NULL, ...)
}
\arguments{
\item{S_list}{List of K sample covariance matrices (each p x p).}
\item{n_vec}{Integer vector of length K with sample sizes.}
\item{burnin}{Integer burn-in iterations.}
\item{nsave}{Integer number of saved posterior draws per chain (after thinning).}
\item{thin}{Integer thinning interval.}
\item{nchains}{Number of independent MCMC chains.}
\item{parallel}{Logical; if TRUE, run chains in parallel using the parallel package.}
\item{ncores}{Number of workers to use when parallel=TRUE.}
\item{seed}{Optional integer seed. If provided, chains use deterministic offsets.}
\item{...}{Additional arguments forwarded to the internal single-chain engine.}
}
\description{
Package-level entry point for the multiGGM MCMC workflow. The current version provides the API and parallel orchestration; the single-chain engine is under active development.
}
\value{
An object of class \code{multiggm_fit} (single chain) or \code{multiggm_fit_list} (multiple chains).
}
