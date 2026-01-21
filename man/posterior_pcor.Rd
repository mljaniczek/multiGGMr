\name{posterior_pcor}
\alias{posterior_pcor}
\title{Extract partial correlations from posterior precision samples}
\usage{
posterior_pcor(Omega_samples)
}
\arguments{
\item{Omega_samples}{Precision samples as a list of matrices or a 3D array (p x p x M).}
}
\description{
Convenience helper to convert posterior samples of precision matrices into partial correlations.
}
\value{
Partial-correlation samples in the same container type as input.
}
