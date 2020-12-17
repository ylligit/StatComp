## ----random, echo=TRUE--------------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
par(mfrow=c(2,2),mar=rep(0,4))
plot(lm.D9)

## ----table, echo=TRUE---------------------------------------------------------
knitr::kable(head(iris))

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------

n <- 10000
a <- 2
b <- 2
unifrv <- runif(n)
paretorv <- b/(1-unifrv)^(1/a)    # inverse transformation
hist(paretorv[paretorv>0 & paretorv<20], freq = FALSE, breaks = seq(0,20,0.5), main = "Histogram of the Pareto sample", xlab = "value")    # graph the density histogram
f <- function(x) {a*b^a/x^(a+1)}    # true pdf
curve(f, 2, 20, col = 2, add = TRUE)    # add the true density curve
legend(12,0.6,"true density", col = 2, lwd = 1)    # add a legend


## -----------------------------------------------------------------------------
u1 <- runif(10000, min = -1, max = 1)    # n = 10000
u2 <- runif(10000, min = -1, max = 1)
u3 <- runif(10000, min = -1, max = 1)
u <- ifelse((abs(u3)>abs(u2) & abs(u3)>abs(u1)), u2, u3)
hist(u, freq = FALSE, breaks = seq(-1,1,0.02), main = "Histogram with the density curve", xlab = "value")
f <- function(x) {3/4*(1-x^2)}
curve(f, -1, 1, col = 2, add = TRUE)   
legend(0.5,0.85,"true density", col = 2, lwd = 1, cex=0.6)    # add a legend


## -----------------------------------------------------------------------------
n <- 1000
r <- 4
beta <- 2
gammarv <- rgamma(n, shape = r, rate = beta)
x <- rexp(n, rate = gammarv)
hist(x[x<5], freq = FALSE, breaks = seq(0,5,0.1), main = "Histogram of the Pareto sample", xlab = "value")  
f <- function(x) {64/(2+x)^5}    # pdf of Pareto distribution
curve(f, 0, 5, col = 2, add = TRUE)    
legend(3, 1, "true density", col = 2, lwd = 1)    # add a legend

## -----------------------------------------------------------------------------
set.seed(1)
m <- 1e5
x <- runif(m, min=0, max=pi/3)
theta.hat <- mean(sin(x)) * pi / 3

## -----------------------------------------------------------------------------
c(theta.hat, 1/2)

## -----------------------------------------------------------------------------
set.seed(2)
m <- 1e6
x1 <- runif(m, min=0, max=1)
theta_hat_1 <- exp(x1)

## -----------------------------------------------------------------------------
set.seed(3)
m <- 5e5
x2 <- runif(m, min=0, max=1);
theta_hat_2 <- (exp(x2) + exp(1-x2)) / 2

## -----------------------------------------------------------------------------
c(1-var(theta_hat_2) / var(theta_hat_1))

## -----------------------------------------------------------------------------
#figure 1
x <- seq(1,5,.01)
g <- x^2/sqrt(2*pi)*exp(-x^2/2)
lambda <- sqrt(2/pi)/exp(1)
f1 <- lambda*exp(-lambda*(x-1))
f2 <- 2/pi/(1+(x-1)^2)
plot(x,g,type="l",main="",ylab="",ylim = c(0,0.4),lwd=2)
lines(x,f1,col='red',lwd=2)
lines(x,f2,col='blue',lwd=2)
legend("topright",legend = c("g","f1","f2"),col=c("black","red","blue"),lwd = 2)

#figure 2
plot(x, g, type = "l", main = "", ylab = "",
ylim = c(0,3), lwd = 2)
lines(x, g/f1, col='red', lwd = 2)
lines(x, g/f2, col='blue', lwd = 2)
legend("topright",legend = c("g","f1","f2"),col=c("black","red","blue"),lwd = 2)

## -----------------------------------------------------------------------------
m <- 1e4
theta.hat <- se <- numeric(2)
g <- function(x){
  x^2/sqrt(2*pi)*exp(-x^2/2)
}
f1 <- function(x,lambda){
  lambda*exp(-lambda*(x-1))
}
f2 <- function(x){
  2/pi/(1+(x-1)^2)
}

#using f1
x <- rexp(m,lambda)+1
fg <- g(x)/f1(x,lambda=lambda)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)
#using f2
u <- runif(m)
x <- tan(pi*u/2)+1
fg <- g(x)/f2(x)
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)

rbind(theta.hat,se)

## -----------------------------------------------------------------------------
u <- seq(0,1,.2)
a <- -log(1-(exp(1)-1)/exp(1)*u)

## -----------------------------------------------------------------------------
inverse.F <- function(u,j){
  -log(1-(exp(1)-1)/(5*exp(1))*(u+j-1))
}
g <- function(x) {
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
f <- function(x){
  exp(-x)/(1-exp(-1))
}

M <- 10000  #number of replicates
k <- 5     #number of strata
r <- M / k  #replicates per stratum
N <- 50     #number of times to repeat the estimation
T2 <- numeric(k)
estimates <- matrix(0,N,2)
set.seed(123)
for (i in 1:N) {
  u <- runif(M)
  x <- - log(1 - u * (1 - exp(-1)))
  fg <- g(x) / f(x)
  estimates[i,1] <- mean(fg)
  for (j in 1:k) {
    u <- runif(r)
    x <- inverse.F(u,j)
    fg <- g(x) / f(x)
    T2[j] <- mean(fg)
  }
  estimates[i,2] <- mean(T2)
}

apply(estimates, 2, mean)
apply(estimates, 2, var)

## -----------------------------------------------------------------------------
n <- 10; m <- 1000
alpha <- .05
CI <- matrix(0,m,2)
set.seed(123)
for (i in 1:m) {
  x <- rnorm(n)
  CI[i,1] <- mean(x)-sd(x)/sqrt(n)*qt(1-alpha/2,n-1)
  CI[i,2] <- mean(x)+sd(x)/sqrt(n)*qt(1-alpha/2,n-1)
}
#confidence level
mean(CI[,1]<0&CI[,2]>0)

## -----------------------------------------------------------------------------
n <- 10; m <- 1000
alpha <- .05
CI <- matrix(0,m,2)
set.seed(123)
for (i in 1:m) {
  x <- rchisq(n,df = 2)
  CI[i,1] <- mean(x)-sd(x)/sqrt(n)*qt(1-alpha/2,n-1)
  CI[i,2] <- mean(x)+sd(x)/sqrt(n)*qt(1-alpha/2,n-1)
}
#confidence level
mean(CI[,1]<2&CI[,2]>2)

## ----echo=FALSE---------------------------------------------------------------

set.seed(12345)

sk = function(x) {
  xbar = mean(x)
  m3 = mean((x - xbar)^3)
  m2 = mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

# beta(a,a)
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

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
 a = c(seq(0,1,0.1),seq(1,20,1),seq(20,100,10))
 pwr = pwr_beta(a)$pwr
 # plot the power
 se = pwr_beta(a)$se
 plot(a, pwr, type = "b", xlab = "a", ylab = "pwr", pch=16)
 abline(h = 0.1, lty = 2)
 lines(a, pwr+se, lty = 4)
 lines(a, pwr-se, lty = 4)

## ----eval=FALSE---------------------------------------------------------------
#  
#  # t(v)
#  pwr_t = function(v){
#  
#   alpha = 0.1
#   n = 20
#   m = 1e3
#   N = length(v)
#   pwr = numeric(N)
#   cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
#  
#   for (j in 1:N) {
#    sktests = numeric(m)
#    for (i in 1:m) {
#     x = rt(n,v[j])
#     sktests[i] = as.integer(abs(sk(x))>= cv)
#    }
#    pwr[j] = mean(sktests)
#   }
#   se = sqrt(pwr*(1-pwr) / m)
#    return(list(pwr = pwr,se = se))
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
v = seq(1,20)
pwr = pwr_t(v)$pwr
se = pwr_t(v)$se
# plot the power
plot(v, pwr, type = "b", xlab = "v", ylab = "pwr", ylim = c(0,1),pch=16)
abline(h = 0.1, lty = 2)
lines(v, pwr+se, lty = 4)
lines(v, pwr-se, lty = 4)


## ----eval=FALSE---------------------------------------------------------------
#  pwr_F <- function(n,alpha.hat){
#    mu1 <- mu2 <- 0
#    sigma1 <- 1
#    sigma2 <- 1.5
#    m <- 1e4
#    result <- matrix(0, length(n), 2)
#    for (i in 1:length(n)){
#      ni <- n[i]
#      tests <- replicate(m, expr={
#        x <- rnorm(ni, mu1, sigma1)
#        y <- rnorm(ni, mu2, sigma2)
#        Fp <- var.test(x, y)$p.value
#        Ftest <- as.integer(Fp <= alpha.hat)
#        c(count5test(x, y), Ftest)
#      })
#      result[i, ] <- rowMeans(tests)
#    }
#    data.frame(n=n, C5=result[, 1], Fp=result[, 2])
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
set.seed(1027)
alpha.hat <- 0.055
n <- c(10, 20, 50, 100, 500, 1000)
pwr <- pwr_F(n,alpha.hat)
print(pwr)

## ----echo=FALSE---------------------------------------------------------------
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

## ----echo=TRUE----------------------------------------------------------------
library(MASS)
library(StatComp20006)
set.seed(1234)
mu <- c(0,0,0)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
m=1000
n<-c(10, 20, 30, 50, 100, 500)
#m: number of replicates; n: sample size
a=numeric(length(n))
for(i in 1:length(n)){
  a[i]=mean(replicate(m, expr={
    mydata <- mvrnorm(n[i],mu,sigma) 
    Mardia(mydata)
  }))
}

## -----------------------------------------------------------------------------
print(a)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(7912)
set.seed(7912)
mu1 <- mu2 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sigma2 <- matrix(c(100,0,0,0,100,0,0,0,100),nrow=3,ncol=3)
sigma=list(sigma1,sigma2)
m=1000
n=50
#m: number of replicates; n: sample size
epsilon <- c(seq(0, .06, .01), seq(.1, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    index=sample(c(1, 2), replace = TRUE, size = n, prob = c(1-e, e))
    mydata<-matrix(0,nrow=n,ncol=3)
    for(t in 1:n){
      if(index[t]==1) mydata[t,]=mvrnorm(1,mu1,sigma1) 
      else mydata[t,]=mvrnorm(1,mu2,sigma2)
    }
    sktests[i] <- Mardia(mydata)
  }
  pwr[j] <- mean(sktests)
}
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## ----eval=FALSE---------------------------------------------------------------
#  jack_correlation <- function(x,y){
#    n <- length(x)
#    theta_hat <- cor(x,y)
#    theta_jack <- numeric(n)
#    for(i in 1:n){
#      theta_jack[i] <- cor(x[-i],y[-i])
#    }
#    bias_jack <- (n-1)*(mean(theta_jack)-theta_hat)
#    se_jack <- sqrt((n-1)*mean((theta_jack-mean(theta_jack))^2))
#    c(bias_jack = bias_jack,se_jack = se_jack)
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
data(law)
attach(law)
jack_est <- jack_correlation(LSAT,GPA)
print(round(jack_est,4))
detach(law)

## ----eval=FALSE---------------------------------------------------------------
#  boot_ci <- function(x){
#    boot.mean <- function(x,i) mean(x[i])
#    boot.obj <- boot(x,statistic = boot.mean,R=2000)
#    boot.ci(boot.obj,type = c("norm","basic","perc","bca"))
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
data(aircondit)
aircondit <- as.matrix(aircondit)
bo <- boot_ci(aircondit)
print(bo)

## ----eval=FALSE---------------------------------------------------------------
#  jack_pca <- function(X){
#    n <- nrow(X)
#    lambda_hat = eigen(cov(X))$values
#    theta_hat = lambda_hat[1] / sum(lambda_hat)
#    theta_j = rep(0,n)
#    for (i in 1:n) {
#      x = X[-i,]
#      lambda = eigen(cov(x))$values
#      theta_j[i] = lambda[1]/sum(lambda)
#    }
#    #estimated bias of theta_hat
#    bias_jack = (n-1)*(mean(theta_j)-theta_hat)
#    #estimated se of theta_hat
#    se_jack = (n-1)*sqrt(var(theta_j)/n)
#    c(bias_jack = bias_jack,se_jack=se_jack)
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
data(scor)
ja <- jack_pca(scor)
print(round(ja,4))

## ----eval=FALSE---------------------------------------------------------------
#  leave2out_cv <- function(X,Y){
#    n <- length(X)
#    e1 <- e2 <- e3 <- e4 <- matrix(0,n*(n-1)/2,2)
#    for(i in 2:n){
#      for(j in 1:(i-1)){
#        k <- (i-1)*(i-2)/2 + j
#        y <- Y[-c(i,j)]
#        x <- X[-c(i,j)]
#  
#        J1 <- lm(y ~ x)
#        yhat1 <- J1$coef[1] + J1$coef[2] * X[c(i,j)]
#        e1[k,] <- Y[c(i,j)]- yhat1
#  
#        J2 <- lm(y ~ x + I(x^2))
#        yhat2 <- J2$coef[1] + J2$coef[2] * X[c(i,j)] + J2$coef[3] * X[c(i,j)]^2
#        e2[k,] <- Y[c(i,j)] - yhat2
#  
#        J3 <- lm(log(y) ~ x)
#        logyhat3 <- J3$coef[1] + J3$coef[2] * X[c(i,j)]
#        yhat3 <- exp(logyhat3)
#        e3[k,] <- Y[c(i,j)] - yhat3
#  
#        J4 <- lm(log(y) ~ log(x))
#        logyhat4 <- J4$coef[1] + J4$coef[2] * log(X[c(i,j)])
#        yhat4 <- exp(logyhat4)
#        e4[k,] <- Y[c(i,j)] - yhat4
#      }
#    }
#    c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
data(ironslag)
attach(ironslag)
cv <- leave2out_cv(chemical,magnetic)
print(cv)
detach(ironslag)

## ----eval=FALSE---------------------------------------------------------------
#  # Count Five test
#  count5test = function(x, y) {
#  X = x - mean(x)
#  Y = y - mean(y)
#  outx = sum(X > max(Y)) + sum(X < min(Y))
#  outy = sum(Y > max(X)) + sum(Y < min(X))
#  # return 1 (reject) or 0 (do not reject H0)
#  return(as.integer(max(c(outx, outy)) > 5))
#  }
#  # Count Five test permutation
#  count5test_permutation <- function(z,R){
#    permutation <- function(z){
#      n <- length(z)
#      x <- z[1:(n/2)]
#      y <- z[-(1:(n/2))]
#      X <- x-mean(x)
#      Y <- y-mean(y)
#      outx <- sum(X > max(Y)) + sum(X < min(Y))
#      outy <- sum(Y > max(X)) + sum(Y < min(X))
#      as.integer(max(c(outx,outy))>5)
#    }
#    n <- length(z)
#    out <- numeric(R)
#    for(r in 1:R){
#      p <- sample(1:n,n,replace = FALSE)
#      out[r] <- permutation(z[p])
#    }
#    sum(out)/R
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
set.seed(1234)
n1 = 20
n2 = 50
mu1 = mu2 = 0
sigma1 = sigma2 = 1
m = 1e3

alphahat1 = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean
y = y - mean(y)
count5test(x, y)
}))
alphahat2 = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean 
y = y - mean(y)
z = c(x,y)
count5test_permutation(z,1000) 
})<0.05)

round(c(count5test=alphahat1,count5test_permutation=alphahat2),4)


## ----eval=FALSE---------------------------------------------------------------
#  eqdist.nn <- function(z,sizes,k,R){
#    Tn <- function(z, ix, sizes,k) {
#      n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#      if(is.vector(z)) z <- data.frame(z,0);
#      z <- z[ix, ];
#      NN <- nn2(data=z, k=k+1) # what's the first column?
#      block1 <- NN$nn.idx[1:n1,-1]
#      block2 <- NN$nn.idx[(n1+1):n,-1]
#      i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#      (i1 + i2) / (k * n)
#    }
#    boot.obj <- boot(data=z,statistic=Tn,R=R,
#                     sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }

## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(Ball)
library(energy)
library(MASS)
library(StatComp20006)

mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0,0,0)
sigma2 <- matrix(c(2,0,0,0,3,0,0,0,4),nrow=3,ncol=3)
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k,R)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,0.5)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k,R)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- as.matrix(rt(n1,1,2),ncol=1)
  mydata2 <- as.matrix(rt(n2,2,5),ncol=1)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k,R)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
rbimodel<-function(n,mu1,mu2,sd1,sd2){
  index=sample(1:2,n,replace=TRUE)
  x=numeric(n)
  index1<-which(index==1)
  x[index1]<-rnorm(length(index1), mu1, sd1)
  index2<-which(index==2)
  x[index2]<-rnorm(length(index2), mu2, sd2)
  return(x)
}
for(i in 1:m){
  mydata1 <- as.matrix(rbimodel(n1,0,0,1,2),ncol=1)
  mydata2 <- as.matrix(rbimodel(n2,1,1,4,3),ncol=1)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k,R)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## -----------------------------------------------------------------------------
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,0.5)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
n1=10
n2=100
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k,R)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

## ----eval=FALSE---------------------------------------------------------------
#  
#  set.seed(3000)
#  
#  lap_f = function(x) exp(-abs(x))
#  
#  rw.Metropolis = function(sigma, x0, N){
#   x = numeric(N)
#   x[1] = x0
#   u = runif(N)
#   k = 0
#   for (i in 2:N) {
#    y = rnorm(1, x[i-1], sigma)
#    if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y
#    else {
#    x[i] = x[i-1]
#    k = k+1
#    }
#   }
#   return(list(x = x, k = k))
#  }
#  
#  N = 2000
#  sigma = c(.05, .5, 2, 16)
#  x0 = 25
#  rw1 = rw.Metropolis(sigma[1],x0,N)
#  rw2 = rw.Metropolis(sigma[2],x0,N)
#  rw3 = rw.Metropolis(sigma[3],x0,N)
#  rw4 = rw.Metropolis(sigma[4],x0,N)
#  #number of candidate points rejected
#  Rej = cbind(rw1$k, rw2$k, rw3$k, rw4$k)
#  Acc = round((N-Rej)/N,4)
#  rownames(Acc) = "Accept rates"
#  colnames(Acc) = paste("sigma",sigma)
#  knitr::kable(Acc)
#  #plot
#  par(mfrow=c(2,2))  #display 4 graphs together
#      rw = cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
#      for (j in 1:4) {
#          plot(rw[,j], type="l",
#               xlab=bquote(sigma == .(round(sigma[j],3))),
#               ylab="X", ylim=range(rw[,j]))
#      }
#  
#  
#  

## ----eval=FALSE---------------------------------------------------------------
#  Gelman.Rubin <- function(psi) {
#  # psi[i,j] is the statistic psi(X[i,1:j])
#  # for chain in i-th row of X
#  psi <- as.matrix(psi)
#  n <- ncol(psi)
#  k <- nrow(psi)
#  psi.means <- rowMeans(psi) #row means
#  B <- n * var(psi.means) #between variance est.
#  psi.w <- apply(psi, 1, "var") #within variances
#  W <- mean(psi.w) #within est.
#  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
#  r.hat <- v.hat / W #G-R statistic
#  return(r.hat)
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
k <- 4    # four chains
x0 <- c(-10,-5,5,10)    # overdispersed initial values
N <- 10000    # length of chains
b <- 200    # burn-in length

opar <- par(mfrow=c(2,2),mar=rep(0,4))
set.seed(123)
X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rwMetropolisR(N,0.5,x0[i])$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (1000+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(1000+1):N], type="l", xlab="sigma=0.5", ylab="R_hat")
abline(h=1.2, lty=2)

X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rwMetropolisR(N,1,x0[i])$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (500+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
x2 <- min(which(rhat>0 & rhat<1.2))
plot(rhat[(500+1):N], type="l", xlab="sigma=1", ylab="R_hat")
abline(h=1.2, lty=2)

X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rwMetropolisR(N,4,x0[i])$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (b+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
x3 <- min(which(rhat>0 & rhat<1.2))
plot(rhat[(b+1):N], type="l", xlab="sigma=4", ylab="R_hat")
abline(h=1.2, lty=2)

X <- matrix(nrow=k,ncol=N)
for (i in 1:k)
  X[i,] <- rwMetropolisR(N,16,x0[i])$x
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
rhat <- rep(0, N)
for (j in (b+1):N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
x4 <- min(which(rhat>0 & rhat<1.2))
plot(rhat[(b+1):N], type="l", xlab="sigma=16", ylab="R_hat")
abline(h=1.2, lty=2)
par(opar)

c(x2,x3,x4)

## ----eval=FALSE---------------------------------------------------------------
#  inter.points <- function(k){
#    S = function(a,k){
#      ck = sqrt(a^2*k/(k+1-a^2))
#      pt(ck,df=k,lower.tail=FALSE)
#    }
#    solve = function(k){
#      output = uniroot(function(a){S(a,k)-S(a,k-1)},lower=1,upper=2)
#      output$root
#    }
#    root = matrix(0,2,length(k))
#    for (i in 1:length(k)){
#      root[2,i]=round(solve(k[i]),4)
#    }
#    root[1,] = k
#    rownames(root) = c('k','A(k)')
#    root
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
k = c(4:25,100,500,1000)
root <- inter.points(k)
print(root)

## ----eval=FALSE---------------------------------------------------------------
#  EM_blood <- function(n.A,n.B,nOO,nAB){
#    # Mle
#    eval_f0 = function(x,x1,n.A,n.B,nOO,nAB) {
#  
#      r1 = 1-sum(x1)
#      nAA = n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
#      nBB = n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
#      r = 1-sum(x)
#      return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
#               (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
#    }
#  
#  
#    # constraint
#    eval_g0 = function(x,x1,n.A,n.B,nOO,nAB) {
#      return(sum(x)-0.999999)
#    }
#  
#    opts = list("algorithm"="NLOPT_LN_COBYLA",
#                "xtol_rel"=1.0e-8)
#    mle = NULL
#    r = matrix(0,1,2)
#    r = rbind(r,c(0.2,0.35))# the beginning value of p0 and q0
#    j = 2
#    while (sum(abs(r[j,]-r[j-1,]))>1e-8) {
#      res = nloptr( x0=c(0.2,0.25),
#                    eval_f=eval_f0,
#                    lb = c(0,0), ub = c(1,1),
#                    eval_g_ineq = eval_g0,
#                    opts = opts, x1=r[j,],n.A=n.A,n.B=n.B,nOO=nOO,nAB=nAB)
#      j = j+1
#      r = rbind(r,res$solution)
#      mle = c(mle,eval_f0(x=r[j,],x1=r[j-1,],n.A=n.A,n.B=n.B,nOO=nOO,nAB=nAB))
#    }
#    return(list(r=r,mle=mle))
#  }

## ----eval=TRUE----------------------------------------------------------------
library(nloptr)
library(StatComp20006)
n.A <- 444; n.B <- 132; nOO <- 361; nAB <- 63
em <- EM_blood(n.A,n.B,nOO,nAB)
#the result of EM algorithm
em$r
#the max likelihood values
plot(-em$mle,type = 'l')

## ----eval=FALSE---------------------------------------------------------------
#  loop_lapply <- function(formulas, data){
#    #1 for loops
#    lo = vector("list", length(formulas))
#    for (i in seq_along(formulas)){
#      lo[[i]] = lm(formulas[[i]], data = data)
#    }
#    la = lapply(formulas, function(x) lm(formula = x, data = data))
#    return(list(lo=lo,la=la))
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
data(mtcars)
attach(mtcars)
formulas = list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
ll <- loop_lapply(formulas,mtcars)
print(ll$lo)
print(ll$la)

## ----eval=FALSE---------------------------------------------------------------
#  sapply_t_p <- function(){
#    trials = replicate(
#      100,
#      t.test(rpois(10, 10), rpois(7, 10)),
#      simplify = FALSE
#    )
#    # anonymous function:
#    p <- sapply(trials, function(x) x[["p.value"]])
#    return(p)
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
set.seed(123)
p <- sapply_t_p()
print(p)

## -----------------------------------------------------------------------------
library(StatComp20006)
data(mtcars); data(faithful)
datalist <- list(mtcars, faithful)
lapply(datalist, function(x) vapply(x, mean, numeric(1)))

## ----eval=FALSE---------------------------------------------------------------
#  mylapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
#    out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
#    if(simplify == TRUE) return(simplify2array(out))
#    unlist(out)
#  }

## -----------------------------------------------------------------------------
mylapply(datalist,mean,numeric(1))

## ----eval= FALSE--------------------------------------------------------------
#  #include <cmath>
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  
#  //[[Rcpp::export]]
#  NumericVector rwMetropolisC(int N, double sigma, double x0) {
#    NumericVector x(N);
#    x[0] = x0;
#    NumericVector u = runif(N);
#    for (int i = 1; i < N;i++ ) {
#      NumericVector y = rnorm(1, x[i-1], sigma);
#      if (u[i] <= exp(abs(x[i-1])-abs(y[0]))){
#        x[i] = y[0];
#      }
#      else {
#        x[i] = x[i-1];
#      }
#    }
#    return(x);
#  }

## ----eval=FALSE---------------------------------------------------------------
#  library(microbenchmark)
#  
#  rwMetropolisR <- function(N,sigma,x0){
#    x <- numeric(N)
#    x[1] <- x0
#    u <- runif(N)
#    k <- 0
#    for (i in 2:N) {
#      y <- rnorm(1,x[i-1],sigma)
#      if(u[i] <= exp(abs(x[i-1])-abs(y))){
#        x[i] <- y
#      }else{
#        x[i] <- x[i-1]
#        k <- k+1
#      }
#    }
#    return(list(x=x,k=k))
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20006)
library(microbenchmark)

x0 = 25; N = 2000; sigma = 2
ts <- microbenchmark(rwR = rwMetropolisR(N,sigma,x0), rwC = rwMetropolisC(N,sigma,x0))
summary(ts)[,c(1,3,5,6)]

## -----------------------------------------------------------------------------
#qqplot
set.seed(123456)
rwR = rwMetropolisR(N,sigma,x0)$x[-(1:500)]
rwC = rwMetropolisC(N,sigma,x0)[-(1:500)]
qqplot(rwR,rwC)
abline(a=0,b=1,col='black')

## ----echo=FALSE---------------------------------------------------------------
knn.est <- function(x, k, xrange, yrange){
p <- ncol(x)
n <- nrow(x)
est_pt <- expand.grid(xrange, yrange)
distance <- knnx.dist(x, est_pt, k)
est_de <- matrix(k / (2 * n * distance[,k]), nrow = length(xrange))
est_de
}

## -----------------------------------------------------------------------------
library(FNN)
library(StatComp20006)
data(faithful)
X <- as.matrix(faithful)
cat(range(X[,1]),range(X[,2]))
xrange <- seq(from=1,to=6,by=.1)
yrange <- seq(from=40,to=100,by=.5)
k <- 5
fit_knnde <- knn.est(X, k, xrange, yrange)
persp(xrange, yrange, fit_knnde, phi = 30, theta = 45, border = 0, col = "blue")

## ----eval=FALSE---------------------------------------------------------------
#  Tstat <- function(x, y){ # return T
#    n <- length(x) # number of samples
#    hx <- hpi(x) # bandwidth for x
#    hy <- hpi(y) # bandwidth for y
#    fx <- rep(0, n)
#    fy <- rep(0, n)
#    # leave-one-out kernel estimators of f(x),f(y)
#    Kx <- matrix(rep(0, n^2),n)
#    Ky <- matrix(rep(0, n^2),n)
#    # The (i,j)-th element of Kx is K_hx(Xj-Xi)
#    # Compute K: Gaussian kernel
#    index <- 1:n
#    for (i in index) {
#      for (j in index[-i]) {
#        Kx[i,j] <- dnorm((x[j] - x[i]) / hx)/ hx
#        Ky[i,j] <- dnorm((y[j] - y[i]) / hy)/ hy
#      }
#    }
#    fx <- n / (n-1) * colMeans(Kx)
#    fy <- n / (n-1) * colMeans(Ky)
#    I <- n / (n-1) * mean(Kx * Ky) + mean(fx %*% t(fy)) - 2 * mean(fx * fy)
#    sigma <- sqrt(2 * sum(Kx^2 * Ky^2) / (n^2 * hx * hy))
#    T <- n * sqrt(hx * hy) * I / sigma
#    return(T)
#  }

## -----------------------------------------------------------------------------
library(StatComp20006)
set.seed(0)
alpha <- 0.05 # alpha for the test
x.true <- as.numeric(faithful[,1]) # x is for eruptions
y.true <- as.numeric(faithful[,2]) # y is for waiting
n <- length(x.true) # number of samples
B <- 20 # bootstrap time

### (1) leave-one-out method
T.true <- Tstat(x.true, y.true)
if (abs(T.true) > qnorm(1 - alpha / 2)){
cat('Using leave-one-out, we reject independent structure.')
} else{
cat('Using leave-one-out, we accept independent structure.')
}

### (2) Boostrap
Tstar <- rep(0,B)
for (i in 1:B) {
xstar <- sample(x.true, n, replace = TRUE)
ystar <- sample(x.true, n, replace = TRUE)
Tstar[i] <- Tstat(xstar, ystar)
} # Bootstrap
p.boot <- mean(T.true < Tstar)
if (p.boot < alpha){
cat('Using bootstrap, we reject independent structure.')
} else{
cat('Using bootstrap, we accept independent structure.')
}

### (3) Compare the two methods
cat('p-value for true T statistics is: ', pnorm(-abs(T.true)) * 2)
cat('p-value for the bootstrap method is: ', p.boot)

