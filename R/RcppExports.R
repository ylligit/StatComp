# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title A random walk Metropolis sampler by using Rcpp
#' @description A random walk Metropolis sampler for generating the standard Laplace distribution by using Rcpp
#' @param N the number of samples
#' @param x0 the initial value
#' @param sigma the standard deviation
#' @return a list contains the generated samplers x of size n and the counter k which records the number of rejected candidate points
#' @examples
#' \dontrun{
#' N <- 2000; sigma <- 2; x0 <- 25
#' rw <- rwMetropolisC(N,sigma,x0)
#' plot(x=1:N,y=rw$x,type="l")
#' }
#' @export
rwMetropolisC <- function(N, sigma, x0) {
    .Call('_StatComp20006_rwMetropolisC', PACKAGE = 'StatComp20006', N, sigma, x0)
}

