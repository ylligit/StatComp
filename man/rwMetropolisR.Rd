% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatCompR.R
\name{rwMetropolisR}
\alias{rwMetropolisR}
\title{A random walk Metropolis sampler by using R}
\usage{
rwMetropolisR(N, sigma, x0)
}
\arguments{
\item{N}{the number of samples}

\item{sigma}{the standard deviation}

\item{x0}{the initial value}
}
\value{
a list contains the generated samplers \code{x} of size \code{n}
and the counter \code{k} which records the number of rejected candidate points.
}
\description{
A random walk Metropolis sampler for generating the standard Laplace distribution by using R
}
\examples{
\dontrun{
N <- 2000; sigma <- 2; x0 <- 25
rw <- rwMetropolisR(N,sigma,x0)
plot(x=1:N,y=rw$x,type="l")
}
}
