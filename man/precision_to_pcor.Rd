\name{precision_to_pcor}
\alias{precision_to_pcor}
\title{Convert precision matrix to partial correlation matrix}
\usage{
precision_to_pcor(Omega)
}
\arguments{
\item{Omega}{A square numeric matrix (precision).}
}
\description{
Given a precision matrix, returns the corresponding partial correlation matrix with diagonal 1.
}
\value{
A numeric matrix of partial correlations.
}
\examples{
Omega <- diag(3)
Omega[1,2] <- Omega[2,1] <- -0.2
precision_to_pcor(Omega)
}
