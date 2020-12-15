#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of R function (\code{rwMetropolisR}) and Cpp function (\code{rwMetropolisC}).
#' @examples
#' \dontrun{
#' x0 = 25; N = 2000; sigma = 2
#' ts <- microbenchmark::microbenchmark(
#'   rwR = rwMetropolisR(N,sigma,x0),
#'   rwC = rwMetropolisC(N,sigma,x0)
#' )
#' print(summary(ts)[,c(1,3,5,6)])
#' }
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @import stats
#' @useDynLib StatComp20006
NULL

#' @title A random walk Metropolis sampler by using R
#' @description A random walk Metropolis sampler for generating the standard Laplace distribution by using R
#' @param N the number of samples
#' @param x0 the initial value
#' @param sigma the standard deviation
#' @return a list contains the generated samplers \code{x} of size \code{n}
#' and the counter \code{k} which records the number of rejected candidate points.
#' @examples
#' \dontrun{
#' N <- 2000; sigma <- 2; x0 <- 25
#' rw <- rwMetropolisR(N,sigma,x0)
#' plot(x=1:N,y=rw$x,type="l")
#' }
#' @export
rwMetropolisR <- function(N,sigma,x0){
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1,x[i-1],sigma)
    if(u[i] <= exp(abs(x[i-1])-abs(y))){
      x[i] <- y
    }else{
      x[i] <- x[i-1]
      k <- k+1
    }
  }
  return(list(x=x,k=k))
}

#' @title Use EM algorithm to solve a MLE problem about A-B-O blood type
#' @description Use EM algorithm to solve MLE of p and q in A-B-O blood type problem condsidering missing data nAA and nBB.
#' @param n.A the number of samples which are A-type.
#' @param n.B the number of samples which are B-type.
#' @param nOO the number of samples which are O-type.
#' @param nAB the number of samples which are AB-type.
#' @return a list contains a matrix r composed by the sequences of p and q value, and a numeric vector MLE which is their correponding log-maximum likelihood.
#' @examples
#' \dontrun{
#' n.A <- 444; n.B <- 132; nOO <- 361; nAB <- 63
#' t <- EM_blood(n.A,n.B,nOO,nAB)
#' t$r
#' plot(-t$mle,type='l')
#' }
#' @importFrom  nloptr nloptr
#' @export
EM_blood <- function(n.A,n.B,nOO,nAB){
  # Mle
  eval_f0 = function(x,x1,n.A,n.B,nOO,nAB) {

    r1 = 1-sum(x1)
    nAA = n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
    nBB = n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
    r = 1-sum(x)
    return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
             (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
  }


  # constraint
  eval_g0 = function(x,x1,n.A,n.B,nOO,nAB) {
    return(sum(x)-0.999999)
  }

  opts = list("algorithm"="NLOPT_LN_COBYLA",
              "xtol_rel"=1.0e-8)
  mle = NULL
  r = matrix(0,1,2)
  r = rbind(r,c(0.2,0.35))# the beginning value of p0 and q0
  j = 2
  while (sum(abs(r[j,]-r[j-1,]))>1e-8) {
    res = nloptr( x0=c(0.2,0.25),
                  eval_f=eval_f0,
                  lb = c(0,0), ub = c(1,1),
                  eval_g_ineq = eval_g0,
                  opts = opts, x1=r[j,],n.A=n.A,n.B=n.B,nOO=nOO,nAB=nAB)
    j = j+1
    r = rbind(r,res$solution)
    mle = c(mle,eval_f0(x=r[j,],x1=r[j-1,],n.A=n.A,n.B=n.B,nOO=nOO,nAB=nAB))
  }
  return(list(r=r,mle=mle))
}

#' @title Use both for loops and \code{lapply()} to fit linear models.
#' @description Use both for loops and \code{lapply()} to fit linear models in data mtcars.
#' @param data from datasets
#' @param formulas a list of models
#' @return A list composed of two lists of fitting models
#' @examples
#' \dontrun{
#' data(mtcars)
#' attach(mtcars)
#' formulas = list(
#'  mpg ~ disp,
#'  mpg ~ I(1 / disp),
#'  mpg ~ disp + wt,
#'  mpg ~ I(1 / disp) + wt
#'  )
#' ll <- loop_lapply(formulas,mtcars)
#' print(ll$lo)
#' print(ll$la)
#' }
#' @export
loop_lapply <- function(formulas, data){
  #1 for loops
  lo = vector("list", length(formulas))
  for (i in seq_along(formulas)){
    lo[[i]] = lm(formulas[[i]], data = data)
  }
  la = lapply(formulas, function(x) lm(formula = x, data = data))
  return(list(lo=lo,la=la))
}

#' @title data
#' @name mtcars
#' @description A dataset used to illustrate the performance of loops and \code{lapply()}.
#' @examples
#' \dontrun{
#' data(mtcars)
#' attach(mtcars)
#' formulas = list(
#'  mpg ~ disp,
#'  mpg ~ I(1 / disp),
#'  mpg ~ disp + wt,
#'  mpg ~ I(1 / disp) + wt
#'  )
#' ll <- loop_lapply(formulas,mtcars)
#' print(ll$lo)
#' print(ll$la)
#' }
NULL

#' @title Use \code{sapply()} to extract the p value from t-test
#' @description The fucntion simulates the performance of a t-test
#'  for non-normal data and use \code{sapply()} and an anonymous
#'   function to extract the p-value from every trial.
#' @return the p values of every trial
#' @examples
#' \dontrun{
#' set.seed(123)
#' p <- sapply_t_p()
#' print(p)
#' }
#' @export
sapply_t_p <- function(){
  trials = replicate(
    100,
    t.test(rpois(10, 10), rpois(7, 10)),
    simplify = FALSE
  )
  # anonymous function:
  p <- sapply(trials, function(x) x[["p.value"]])
  return(p)
}

#' @title A combination of \code{Map()} and \code{vapply()} to create an \code{lapply()} variant.
#' @description Implement a combination of \code{Map()} and \code{vapply()}
#' to create an \code{lapply()} variant that iterates in parallel over all of its inputs
#' and stores its outputs in a vector (or a matrix).
#' @param X a vector (atomic or list) or an expression object.
#' @param FUN the function to be applied to each element of \code{X}
#' @param FUN.VALUE a (generalized) vector; a template for the return value from FUN.
#' @param simplify ogical or character string; should the result be simplified to a vector,
#' matrix or higher dimensional array if possible.
#' @return a vector or a matrix of the value generated by applying FUN to the given inputs \code{X}
#' @examples
#' \dontrun{
#'  data(mtcars); data(faithful)
#'  datalist <- list(mtcars,faithful)
#'  mylapply(datalist, mean, numeric(1))
#' }
#' @export
mylapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE) return(simplify2array(out))
  unlist(out)
}

#' @title data
#' @name faithful
#' @description a dataset from package datasets
NULL

#' @title Use the Gelman-Rubin method to monitor convergence of the chain generated by \code{rwMetropolisR}.
#' @description Use the Gelman-Rubin method to monitor convergence of the chain
#'  generated by \code{rwMetropolisR}, and run the chain until it converges
#'  approximately to the target distribution according to \code{R_hat <1.2} .
#' @param psi \code{psi[i,j]} is the statistic \code{psi(X[i,1:j])}
#' @return a numeric value of G-R statistic
#' @examples
#' \dontrun{
#'  library(StatComp20006)
#'  k <- 4    # four chains
#'  x0 <- c(-10,-5,5,10)    # overdispersed initial values
#'  N <- 10000    # length of chains
#'  b <- 200    # burn-in length
#'  X <- matrix(nrow=k,ncol=N)
#'  for (i in 1:k)
#'    X[i,] <- rwMetropolisR(N,1,x0[i])$x
#'  psi <- t(apply(X, 1, cumsum))
#'  for (i in 1:nrow(psi))
#'    psi[i,] <- psi[i,] / (1:ncol(psi))
#'  rhat <- rep(0, N)
#'  for (j in (500+1):N)
#'    rhat[j] <- Gelman.Rubin(psi[,1:j])
#'  x2 <- min(which(rhat>0 & rhat<1.2))
#'  plot(rhat[(500+1):N], type="l", xlab="sigma=1", ylab="R_hat")
#'  abline(h=1.2, lty=2)
#' }
#' @export
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

#' @title Find the intersection points of a kind of curves
#' @description Find the intersection points of A(k) in (0,sqrt(k)) of the curves
#' S_{k-1}(a) = P(t(k-1)>sqrt{frac{a^2(k-1)}{k-a^2}}) and
#' S_k(a) = P(t(k)>sqrt{frac{a^2k}{k+1-a^2}}),
#' where t(k) is a Student t random variable with k degrees of freedom.
#' @param k A numeric vector
#' @return A matrix of size 2*length(k). The first row is 1:k
#' and the second row is the corresponding root value.
#' @examples
#' \dontrun{
#' k = c(4:25,100,500,1000)
#' root <- inter.points(k)
#' print(root)
#' }
#' @export
inter.points <- function(k){
  S = function(a,k){
    ck = sqrt(a^2*k/(k+1-a^2))
    pt(ck,df=k,lower.tail=FALSE)
  }
  solve = function(k){
    output = uniroot(function(a){S(a,k)-S(a,k-1)},lower=1,upper=2)
    output$root
  }
  root = matrix(0,2,length(k))
  for (i in 1:length(k)){
    root[2,i]=round(solve(k[i]),4)
  }
  root[1,] = k
  rownames(root) = c('k','A(k)')
  root
}

#' @title The two sample “Count Five” test for equality of variance.
#' @description Suppose the means of the two samples are equal and the sample sizes are equal.
#'  An observation in one sample is considered extreme if it is not within the range of the other sample.
#'  If either sample has five or more extreme points, the hypothesis of equal variance is rejected.
#' @param x a numeric vector which represents a sample
#' @param y a numeric vector which represents another sample
#' @return 1 represents that the hypothesis is rejected, while 0 represents that the hypothesis is not rejected.
#' @examples
#' \dontrun{
#' n1 <- n2 <- 20
#' mu1 <- mu2 <- 0
#' sigma1 <- 1; sigma2 <- 2
#' m <- 1e3
#' alphahat <- mean(replicate(m,expr={
#' x=rnorm(n1,mu1,sigma1)
#' y=rnorm(n2,mu2,sigma2)
#' count5test(x,y)
#' }))
#' print(round(alphahat,4))
#' }
#' @export
count5test = function(x, y) {
  X = x - mean(x)
  Y = y - mean(y)
  outx = sum(X > max(Y)) + sum(X < min(Y))
  outy = sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

#' @title Count five permutation test
#' @description Implement a permutation test for equal variance
#' based on the maximum number of extreme points that applies
#' when sample sizes are not necessarily equal.
#' @param z a numeric vector which is generated by combining sample x and sample y
#' @param R the permutation times
#' @return a number in [0,1] which represents the frequency of count5test being rejected.
#' @examples
#' \dontrun{
#' n1 <- 20; n2 <-50
#' mu1 <- mu2 <- 0
#' sigma1 <- sigma2 <- 1
#' m <- 1e3
#' alphahat <- mean(replicate(m,expr={
#' x=rnorm(n1,mu1,sigma1)
#' y=rnorm(n2,mu2,sigma2)
#' z=c(x,y)
#' count5test_permutation(z,1000)}))
#' print(alphahat)
#' }
#' @export
count5test_permutation <- function(z,R){
  permutation <- function(z){
    n <- length(z)
    x <- z[1:(n/2)]
    y <- z[-(1:(n/2))]
    X <- x-mean(x)
    Y <- y-mean(y)
    outx <- sum(X > max(Y)) + sum(X < min(Y))
    outy <- sum(Y > max(X)) + sum(Y < min(X))
    as.integer(max(c(outx,outy))>5)
  }
  n <- length(z)
  out <- numeric(R)
  for(r in 1:R){
    p <- sample(1:n,n,replace = FALSE)
    out[r] <- permutation(z[p])
  }
  sum(out)/R
}

#' @title Nearest neighbor test for equal distribution.
#' @description Nearest neighbor test for equal distribution.
#' @param z a numeric vector which represents data sample.
#' @param sizes a numeric vector of length 2
#' @param k a scalar
#' @param R the number of bootstraps
#' @return a list contains the bootstrap statistic and its p value
#' @import boot
#' @import RANN
#' @examples
#' \dontrun{
#' m <- 1e2; k <- 3; p <- 2; mu <- .6
#' n1 <- n2 <- 20; R <- 999; N <- c(n1,n2)
#' p.value <- numeric(m)
#' for(i in 1:m){
#' x <- matrix(rnorm(n1*p),ncol=p)
#' y <- matrix(rnorm(n2*p,mean=mu),ncol=p)
#' z <- rbind(x,y)
#' p.values[i] <- eqdist.nn(z,sizes=N,k=k,R=R)$p.value
#' }
#' alpha <- 0.1
#' pow <- colMeans(p.values<alpha)
#' pow
#' }
#' @export
eqdist.nn <- function(z,sizes,k,R){
  Tn <- function(z, ix, sizes,k) {
    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
    if(is.vector(z)) z <- data.frame(z,0);
    z <- z[ix, ];
    NN <- nn2(data=z, k=k+1) # what's the first column?
    block1 <- NN$nn.idx[1:n1,-1]
    block2 <- NN$nn.idx[(n1+1):n,-1]
    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
    (i1 + i2) / (k * n)
  }
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

#' @title Compute a jackknife estimate of the bias and the standard error
#' of the correlation statistic
#' @description Compute a jackknife estimate of the bias and the standard error
#' of the correlation statistic
#' @param x a numeric vector of sample
#' @param y a numeric vector of another sample
#' @return a vector of size 1*2 which contains the jackknife estimate of bias and
#' the jackknife estimate of the standard error of the correlation statistic.
#' @examples
#' \dontrun{
#' data(law)
#' attach(law)
#' jack_est <- jack_correlation(LSAT,GPA)
#' round(jack_est,4)
#' }
#' @export
jack_correlation <- function(x,y){
  n <- length(x)
  theta_hat <- cor(x,y)
  theta_jack <- numeric(n)
  for(i in 1:n){
    theta_jack[i] <- cor(x[-i],y[-i])
  }
  bias_jack <- (n-1)*(mean(theta_jack)-theta_hat)
  se_jack <- sqrt((n-1)*mean((theta_jack-mean(theta_jack))^2))
  c(bias_jack = bias_jack,se_jack = se_jack)
}

#' @title data
#' @name law
#' @description a dataset from bootstrap package
NULL

#' @title Compute bootstrap confidence intervals for the mean by the standard normal, basic, percentile, and BCa methods.
#' @description Compute 95 bootstrap confidence intervals for the mean by the standard normal, basic, percentile, and BCa methods.
#' @param x a numeric vector of data sample
#' @return a list contains the 4 C.I. and other information.
#' @import boot
#' @examples
#' \dontrun{
#' data(aircondit)
#' aircondit <- as.matrix(aircondit)
#' bo <- boot_ci(aircondit)
#' print(bo)
#' }
#' @export
boot_ci <- function(x){
  boot.mean <- function(x,i) mean(x[i])
  boot.obj <- boot(x,statistic = boot.mean,R=2000)
  boot.ci(boot.obj,type = c("norm","basic","perc","bca"))
}

#' @title data
#' @name aircondit
#' @description a dataset from boot package
NULL

#' @title data
#' @name scor
#' @description a data matrix from bootstrap package
NULL

#' @title Obtain the jackknife estimates of bias and standard
#' error of the proportion of variance explained by the first principal component.
#' @description Obtain the jackknife estimates of bias and standard
#' error of the proportion of variance explained by the first principal component.
#' @param X a data matrix
#' @return a numeric vector contains the jackknife estimate of bias and se.
#' @examples
#' \dontrun{
#' data(scor)
#' ja <- jack_pca(scor)
#' print(round(ja,4))
#' }
#' @export
jack_pca <- function(X){
  n <- nrow(X)
  lambda_hat = eigen(cov(X))$values
  theta_hat = lambda_hat[1] / sum(lambda_hat)
  theta_j = rep(0,n)
  for (i in 1:n) {
    x = X[-i,]
    lambda = eigen(cov(x))$values
    theta_j[i] = lambda[1]/sum(lambda)
  }
  #estimated bias of theta_hat
  bias_jack = (n-1)*(mean(theta_j)-theta_hat)
  #estimated se of theta_hat
  se_jack = (n-1)*sqrt(var(theta_j)/n)
  c(bias_jack = bias_jack,se_jack=se_jack)
}

#' @title leave_two_out cross validation
#' @description Use leave-two-out cross validation to compare the models.
#' @param X a data sample which represents the response varirable.
#' @param Y another data sample which represents the predictive variable.
#' @return a numeric vector of residual sum of squares of all models.
#' @examples
#' \dontrun{
#' data(ironslag)
#' attach(ironslag)
#' cv <- leave2out_cv(chemical,magnetic)
#' print(cv)
#' }
#' @export
leave2out_cv <- function(X,Y){
  n <- length(X)
  e1 <- e2 <- e3 <- e4 <- matrix(0,n*(n-1)/2,2)
  for(i in 2:n){
    for(j in 1:(i-1)){
      k <- (i-1)*(i-2)/2 + j
      y <- Y[-c(i,j)]
      x <- X[-c(i,j)]

      J1 <- lm(y ~ x)
      yhat1 <- J1$coef[1] + J1$coef[2] * X[c(i,j)]
      e1[k,] <- Y[c(i,j)]- yhat1

      J2 <- lm(y ~ x + I(x^2))
      yhat2 <- J2$coef[1] + J2$coef[2] * X[c(i,j)] + J2$coef[3] * X[c(i,j)]^2
      e2[k,] <- Y[c(i,j)] - yhat2

      J3 <- lm(log(y) ~ x)
      logyhat3 <- J3$coef[1] + J3$coef[2] * X[c(i,j)]
      yhat3 <- exp(logyhat3)
      e3[k,] <- Y[c(i,j)] - yhat3

      J4 <- lm(log(y) ~ log(x))
      logyhat4 <- J4$coef[1] + J4$coef[2] * log(X[c(i,j)])
      yhat4 <- exp(logyhat4)
      e4[k,] <- Y[c(i,j)] - yhat4
    }
  }
  c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
}

#' @title data
#' @name ironslag
#' @description a data matrix from DAAG package
NULL

#' @title Calculate the sample skewness
#' @description Calculate the sample skewness
#' @param x a numeric vector of sample
#' @return a scalar of sample skewness
#' @export
sk = function(x) {
  xbar = mean(x)
  m3 = mean((x - xbar)^3)
  m2 = mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

#' @title Estimate the power of the skewness test of
#' normality against symmetric Beta(a,a) distributions
#' @description Estimate the power of the skewness test of
#' normality against symmetric Beta(a,a) distributions
#' @param a the parametric of Beta distribution
#' @return a list contains the power and its se
#' @examples
#' \dontrun{
#' a <- c(seq(0,1,0.1),seq(1,20,1),seq(20,100,10))
#' pwr <- pwr_beta(a)$pwr
#' se <- pwr_beta(a)$se
#' plot(a,pwr,type="b",xlab="a",ylab="pwr",pch=16)
#' abline(h=0.1,lty=2)
#' lines(a, pwr+se, lty = 4)
#' lines(a, pwr-se, lty = 4)
#' }
#' @export
pwr_beta = function(a){
  alpha = 0.1
  n = 20
  m = 1e4
  N = length(a)
  pwr = numeric(N)
  cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

  for (j in 1:N) {
    sktests = numeric(m)
    for (i in 1:m) {
      x = rbeta(n, a[j], a[j])
      sktests[i] = as.integer(abs(sk(x))>= cv)
    }
    pwr[j] = mean(sktests)
  }
  se = sqrt(pwr * (1-pwr) / m)
  return(list(pwr = pwr,se = se))
}


#' @title Estimate the power of the skewness test of
#' normality against t(v) distributions
#' @description Estimate the power of the skewness test of
#' normality against t(v) distributions
#' @param v the parametric of t(v) distributions.
#' @return a list contains the power and its se
#' @examples
#' \dontrun{
#' v <- seq(1,20)
#' pwr <- pwr_t(v)$pwr
#' se <- pwr_t(v)$se
#' plot(v,pwr,type="b",xlab="v",ylab="pwr",ylim=c(0,1),pch=16)
#' abline(h=0.1,lty=2)
#' lines(v, pwr+se, lty = 4)
#' lines(v, pwr-se, lty = 4)
#' }
#' @export
pwr_t = function(v){

  alpha = 0.1
  n = 20
  m = 1e3
  N = length(v)
  pwr = numeric(N)
  cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

  for (j in 1:N) {
    sktests = numeric(m)
    for (i in 1:m) {
      x = rt(n,v[j])
      sktests[i] = as.integer(abs(sk(x))>= cv)
    }
    pwr[j] = mean(sktests)
  }
  se = sqrt(pwr*(1-pwr) / m)
  return(list(pwr = pwr,se = se))
}

#' @title Compute the F test of equal variance.
#' @description Compute the F test of equal variance at significance level 0.055 and compare the power of the Count Five test and F test for small, medium, and large sample sizes.
#' @return a dataframe contains n and power of Count 5 test and F test.
#' @param n a numeric vector
#' @param alpha.hat the confidence level
#' @examples
#' \dontrun{
#' set.seed(1027)
#' alpha.hat <- 0.055
#' n <- c(10, 20, 50, 100, 500, 1000)
#' pwr <- pwr_F(n,alpha.hat)
#' print(pwr)
#' }
#' @export
pwr_F <- function(n,alpha.hat){
  mu1 <- mu2 <- 0
  sigma1 <- 1
  sigma2 <- 1.5
  m <- 1e4
  result <- matrix(0, length(n), 2)
  for (i in 1:length(n)){
    ni <- n[i]
    tests <- replicate(m, expr={
      x <- rnorm(ni, mu1, sigma1)
      y <- rnorm(ni, mu2, sigma2)
      Fp <- var.test(x, y)$p.value
      Ftest <- as.integer(Fp <= alpha.hat)
      c(count5test(x, y), Ftest)
    })
    result[i, ] <- rowMeans(tests)
  }
  data.frame(n=n, C5=result[, 1], Fp=result[, 2])
}

#' @title Mardia's multivariate skewness test.
#' @description Mardia's multivariate skewness test.
#' @param mydata a numeric matrix of data
#' @return 0 or 1, 1 represents the test being rejected, which 0 represents the test cannot be rejected.
#' @examples
#' \dontrun{
#' mu <- c(0,0,0)
#' sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
#' m=1000
#' n<-c(10, 20, 30, 50, 100, 500)
#' a=numeric(length(n))
#' for(i in 1:length(n)){
#'   a[i]=mean(replicate(m, expr={
#'     mydata <- mvrnorm(n[i],mu,sigma)
#'     Mardia(mydata)
#'   }))
#' }
#' }
#' @export
Mardia<-function(mydata){
  n=nrow(mydata)
  c=ncol(mydata)
  central<-mydata
  for(i in 1:c){
    central[,i]<-mydata[,i]-mean(mydata[,i])
  }
  sigmah<-t(central)%*%central/n
  a<-central%*%solve(sigmah)%*%t(central)
  b<-sum(colSums(a^{3}))/(n*n)
  test<-n*b/6
  chi<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>chi)
}

#' @title Use function knn.dist() in package FNN
#' to write a density estimate function
#' @description Use function knn.dist() in package FNN
#' to write a density estimate function
#' @param x an input data matrix
#' @param k the maximum number of nearest neighbors to search.
#' @param xrange the range of x sample
#' @param yrange the range of y sample
#' @return a matrix contains the value of estimation
#' @importFrom FNN knnx.dist
#' @examples
#' \dontrun{
#' data(faithful)
#' X <- as.matrix(faithful)
#' xrange <- seq(from=1,to=6,by=.1)
#' yrange <- seq(from=40,to=100,by=.5)
#' k <- 5
#' fit_knnde <- knn.est(X, k, xrange, yrange)
#' persp(xrange, yrange, fit_knnde, phi = 30, theta = 45, border = 0, col = "blue")
#' }
#' @export
knn.est <- function(x, k, xrange, yrange){
  p <- ncol(x)
  n <- nrow(x)
  est_pt <- expand.grid(xrange, yrange)
  distance <- knnx.dist(x, est_pt, k)
  est_de <- matrix(k / (2 * n * distance[,k]), nrow = length(xrange))
  est_de
}

#' @title Generate random data from a mixnormal distribution
#' @description Generate n data points from the distribution
#' 0.3 N(0, 1) + 0.7 N(1, 0.32), use the bandwidth selection methods
#' in R package kedd, and draw the estimated density curves with
#' different bandwidths in one plot.
#' @param n a scalar which represents the amount of generated data points.
#' @return the generated data
#' @import kedd
#' @importFrom graphics legend lines
#' @examples
#' \dontrun{
#' n <- 100
#' set.seed(0)
#' mixnorm(n)
#' }
#' @export
mixnorm <- function(n){
  dMG <- function(x) 0.3 * dnorm(x,0,1) + 0.7 * dnorm(x,1,0.3)
  rMG <- function(n){
    # randomly generate n points from the Mixed Gaussian distribution
    r <- runif(n, 0, 1)
    x <- r
    ind <- which(r < 0.3) #index for those generated from N(0,1)
    x[ind] <- rnorm(length(ind), 0, 1)
    x[-ind] <- rnorm(n-length(ind), 1, 0.3)
    return(x)
  }
  x <- rMG(n)
  fhat.amise <- dkde(x, h = h.amise(x)$h) # BW selection: amise
  fhat.bcv <- dkde(x, h = h.bcv(x)$h) # BW selection: bcv
  fhat.ccv <- dkde(x, h = h.ccv(x)$h) # BW selection: ccv

  plot(dMG, from = -2, to = 3, lwd = 2)
  lines(fhat.amise$eval.points, fhat.amise$est.fx, col = 2, lwd = 2)
  lines(fhat.bcv$eval.points, fhat.bcv$est.fx, col = 3, lwd = 2)
  lines(fhat.ccv$eval.points, fhat.ccv$est.fx, col = 4, lwd = 2)
  legend('topleft', legend = c('true', 'amise', 'bcv', 'ccv'),
         col = 1:4, lwd = c(2, 2, 2, 2))
  return(x)
}
