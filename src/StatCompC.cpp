#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk Metropolis sampler by using Rcpp
//' @description A random walk Metropolis sampler for generating the standard Laplace distribution by using Rcpp
//' @param N the number of samples
//' @param x0 the initial value
//' @param sigma the standard deviation
//' @return a list contains the generated samplers x of size n and the counter k which records the number of rejected candidate points
//' @examples
//' \dontrun{
//' N <- 2000; sigma <- 2; x0 <- 25
//' rw <- rwMetropolisC(N,sigma,x0)
//' plot(x=1:N,y=rw$x,type="l")
//' }
//' @export
// [[Rcpp::export]]
NumericVector rwMetropolisC(int N, double sigma, double x0) {
  NumericVector x(N);
  x[0] = x0;
  NumericVector u = runif(N);
  for (int i = 1; i < N;i++ ) {
    NumericVector y = rnorm(1, x[i-1], sigma);
    if (u[i] <= exp(abs(x[i-1])-abs(y[0]))){
      x[i] = y[0];
    }
    else {
      x[i] = x[i-1];
    }
  }
  return(x);
}
