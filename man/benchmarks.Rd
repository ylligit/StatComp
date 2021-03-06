% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatCompR.R
\name{benchmarks}
\alias{benchmarks}
\title{Benchmark R and Rcpp functions.}
\description{
Use R package \code{microbenchmark} to compare the performance of R function (\code{rwMetropolisR}) and Cpp function (\code{rwMetropolisC}).
}
\examples{
\dontrun{
x0 = 25; N = 2000; sigma = 2
ts <- microbenchmark::microbenchmark(
  rwR = rwMetropolisR(N,sigma,x0),
  rwC = rwMetropolisC(N,sigma,x0)
)
print(summary(ts)[,c(1,3,5,6)])
}
}
