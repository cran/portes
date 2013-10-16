\name{portes-package}
\alias{portes-package}
\docType{package}

\title{
Portmanteau Tests for Univariate and Multivariate Time Series Models 
}

\description{
This package contains a set of portmanteau diagnostic checks for univariate and multivariate time series
based on the asymptotic approximation distributions and the Monte-Carlo significance test.
More details about the portmanteau test statistics are given in the online vignette of this package.
It can be used for generating a simulated data from nonseasonal \code{ARIMA} or \code{VARIMA} models whith
innovations from finite or infinite variances distributions. 
The simulated data may have deterministic terms, a constant drift and a time trend, with non-zero mean.
}

\details{
\tabular{ll}{
Package: \tab portes\cr
Type: \tab Package\cr
Version: \tab 2.1-1\cr
Date: \tab 2013-10-16\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
Depends: \tab \code{R (>= 2.14.0), parallel}\cr 
Suggests: \tab \code{fGarch(V.2150.81), FitAR(V.1.92), FGN(V.1.5), TSA(V.0.99), }\cr
\tab \code{vars(V.1.5-0), tseries(V.0.10-29), forecast(V.3.25), akima(V.0.5-7)}\cr
Classification/ACM: \tab G.3, G.4, I.5.1\cr
Classification/MSC: \tab 62M10, 91B84\cr
License: \tab GPL (>= 2)\cr
}
}
\section{Main Function}{
The main function in this package, \code{\link{portest}}, 
is used with univariate and multivariate time series.
It implements the Portmanteau test statistics, \code{\link{gvtest}}, 
\code{\link{BoxPierce}}, \code{\link{LjungBox}}, \code{\link{Hosking}}, 
and \code{\link{LiMcLeod}} based on two methods.
The first method uses the Monte-Carlo techniques as described by 
Lin and McLeod (2006), Mahdi and McLeod (2012) and the second one uses 
the approximation asymptotic distribution.
Originally, the generalized variance portmanteau test, \code{\link{gvtest}}, 
for univariate time series was derived by Pena and Rodriguez (2002)
based on the gamma distribution. 
Lin and McLeod (2006) proposed the Monte-Carlo version of this test and 
Mahdi and McLeod (2012) extended both methods to the multivariate case. 
Simulation results suggest that the Monte-Carlo version of 
\code{\link{gvtest}} statistic is more accurate
and powerful than its competitors proposed by Box and Pierce (1970), Ljung and Box (1978),
and Pena and Rodriguez (2002, 2006) in the univariate time series and
Hosking (1980) and Li and McLeod (1981) in the multivariate time series.

The powerful parallel computing framework facility
is implemented in this function using the package \code{parallel}.
This package handles running much larger chunks of computations in parallel and was first 
included in \code{R 2.14.0} based on the work done for \code{CRAN} packages
\code{multicore} and \code{snow}.

The default argument in \code{\link{portest}} function, \code{nslaves=1}, 
implements the Monte-Carlo test on PC with only one CPU. 
Set the argument \code{nslaves} equals to a positive integer number greater than 1, 
provided that the default argument \code{MonteCarlo=TRUE} is selected, 
the package \code{parallel} requires the Message Passing Interface, MPI, language 
to be properly installed and implemented by the parallel computing program \code{MPICH2}    
in the I O environment system in which the \code{portes} package will run.
Instructions to install and run \code{MPICH2} is given in the link
\url{http://www.stats.uwo.ca/faculty/yu/Rmpi}.
When \code{MonteCarlo=FALSE} is selected, the test statistic 
selected from the argument \code{test} will be implemented
based on the asymptotic approximation distribution.
The default test statistic is the generalized variance test, \code{\link{gvtest}}.

\subsection{Test for usual residuals and \code{GARCH} effects}{

By setting the argument \code{SquaredQ=TRUE} in \code{\link{portest}} function,
the portmanteau test using the asymptotic distribution approximation or the 
Monte-Carlo significance test 
(depending on the choice of the argument \code{MonteCarlo} whether \code{FALSE} or \code{TRUE})
for \code{ARCH} effects will be implemented on the squared residuals.
Otherwise, the portmanteau test will be applied on the usual residuals (when the 
default argument \code{SquaredQ=FALSE} is selected).
The goodness-of-fit garch model test can be implemented on fitted models 
with classes \code{"garch"} and \code{"fGarch"}. 
These two classes are associated with output objects from the
functions \code{garch()} and \code{garchFit()}
available from the R packages \code{tseries} and \code{fGarch} respectively.
}

\subsection{Monte-Carlo test for residuals with infinite variances}{

The argument \code{InfiniteVarianceQ=TRUE} in \code{\link{portest}} function is used only with
Monte-Carlo techniques. 
By selecting this argument, the Monte-Carlo diagnostic test on residuals with
infinite variances effects is implemented.
}

\subsection{Test for fractional Gaussian noise, \code{FGN}, effects}{

After fitting \code{FGN} model using the function \code{FitFGN()} available from the \code{FGN} R package,
the output object has a class \code{"FitFGN"}.
By substituting this object as a first entry in the \code{\link{portest}} function, the portmanteau
test based on the Monte-Carlo or the asymptotic distribution method 
(depending on the choice of the argument \code{MonteCarlo} whether \code{FALSE} or \code{TRUE})
for \code{FGN} model will be implemented.
}

\subsection{Goodness-of-fit test for any fitted model}{
The portmanteau test statistics implemented in the functions, \code{\link{gvtest}}, \code{\link{BoxPierce}}, 
\code{\link{LjungBox}}, \code{\link{Hosking}}, \code{\link{LiMcLeod}}, can be 
used for testing the adequacy of any fitted model using the Monte-Carlo significance test 
or the asymptotic aprroximation distribution test. 
More details with illustrative examples are given in the documentation of
the main function \code{\link{portest}}.
}
}

\section{Simulate data from nonseasonal \code{ARIMA(p,d,q)} or \code{VARIMA(p,d,q)}}{
The function \code{\link{varima.sim}} in this package is useful for simulating data
from nonseasonal \code{ARIMA} or \code{VARIMA} of order \eqn{(p,d,q)} 
with or without deterministic terms (drift and trend).
The innovations series can be given via the argument \code{innov} or
may be generated from the selected distribution specified from the argument \code{innov.dist}.
The argument \code{d} must entered as a nonnegative integer in the \code{ARIMA} case, whereas 
it must entered as a vector of \eqn{k} components \eqn{d_1,...,d_k} in the \code{VARIMA} case.
\eqn{d_i} represents the difference lag need to be applied on series \eqn{i}.
The components of the argument \code{StableParameters} are the 
stable parameters \code{ALPHA, BETA, GAMMA}, and \code{DELTA} needed to generate  
innovations from stable distribution. 
}

\examples{
\dontrun{
#########################################################################
####                                                                 ####
####                   Monte-Carlo Portmanteau Tests                 ####  
####                                                                 #### 
#########################################################################
## Monte-Carlo test for randomness series                              ##
#########################################################################
data("DEXCAUS")
returns <- log(DEXCAUS[-1]/DEXCAUS[-length(DEXCAUS)])
portest(returns)            ## MC using one CPU takes about 25.16 seconds
portest(returns, nslaves=4) ## MC using 4 CPUs takes about 9.51 seconds
portest(returns, MonteCarlo=FALSE)                 ## asymptotic gvtest
portest(returns,test="LjungBox", MonteCarlo=FALSE) ## asymptotic LjungBox 
#########################################################################
## Monte-Carlo goodness-of-fit arima test using 4 CPUs                 ##
#########################################################################
## arima() or Arima() function takes about 14.32 seconds
ans1 <- arima(WWWusage,order=c(3,1,0))
portest(ans1, nslaves = 4)
#
## arima0() function takes about 15.26 seconds
ans2 <- arima0(WWWusage,order=c(3,1,0))
portest(ans2, nslaves = 4)
#
## auto.arima() function from package forecast takes about 13.59 seconds
library("forecast")
ans3 <- auto.arima(WWWusage)
portest(ans3, nslaves = 4)
#
## ar() function takes about 9.39 seconds
ans4 <- ar(Nile,order.max=2)
portest(ans4, nslaves = 4)
#
## FitAR() function takes about 10.78 seconds
library("FitAR")
ans5 <- FitAR(Nile, p=2)
portest(ans5, nslaves = 4)
#########################################################################
## Monte-Carlo goodness-of-fit VAR test - Multivariate series     ##
#########################################################################
data("IbmSp500")
ibm <- log(IbmSp500[,2]+1)*100
sp500 <- log(IbmSp500[,3]+1)*100
IBMSP500 <- data.frame(cbind(ibm,sp500))
## ar.ols() function takes about 9.11 seconds
ans6 <-  ar.ols(IBMSP500, aic=FALSE, intercept=TRUE, order.max=5)
portest(ans6, NREP=100, test="gvtest", nslaves=4)
## VAR() function takes about 11.55 seconds
library("vars")
ans7 <- VAR(IBMSP500, p=5) 
portest(ans7, NREP=100, test="gvtest", nslaves=4)
portest(ans7,test="Hosking", MonteCarlo=FALSE) ## asymptotic Hosking test
#########################################################################
## Monte-Carlo test for GARCH effects using 4 CPUs                     ##
#########################################################################
## Example 1
## Test for ARCH effects on returns series takes about 14.65 seconds
data("monthintel")
returns <- as.ts(monthintel)
lags <- c(5, 10, 20, 40)
portest(returns, lags = lags, nslaves = 4, SquaredQ = TRUE)
#
## Example 2
library("fGarch")
library("tseries")
data("GNPDEF")
z<-ts(GNPDEF[,2], start=1947, freq=4)
r <- 100*diff(log(z))
## use garch() function takes about 6.75 seconds
FitGarch1 <- garch(r, order = c(1,1))
portest(FitGarch1,NREP=100,nslaves = 4)
portest(FitGarch1,NREP=100,nslaves = 4,SquaredQ=TRUE)
#
## use garchFit() function takes about 13.56 seconds
GarchFit2 <- garchFit(formula = ~arma(4,0)+garch(1,1), data=r, trace=FALSE)
portest(GarchFit2, NREP=100, nslaves = 4, SquaredQ = FALSE)
portest(GarchFit2, NREP=100, nslaves = 4, SquaredQ = TRUE)
#########################################################################
## Monte-Carlo test on residuals with infinite variances               ##
#########################################################################
## It takes about 32.7 seconds on personal PC with 4 CPUs 
data("CRSP")
CRSP.AR5<- arima(CRSP, c(5, 0, 0))
lags <- c(10, 20, 30)
portest(CRSP.AR5,lags=lags,nslaves=4,NREP=1000,InfiniteVarianceQ=TRUE)      
#########################################################################
## Monte-Carlo test for Fractional Gaussian Noise, FGN.                ##
#########################################################################
## It takes about 55.06 seconds on personal PC with 4 CPUs 
library("FGN")
data("NileMin")
NILE.FGN <- FitFGN(NileMin)
lags <- c(5, 10, 20)
portest(NILE.FGN, lags = lags, nslaves = 4) 
portest(NILE.FGN, MonteCarlo=FALSE)  ## asymptotic distribution method
##############################################################
## Write two functions to fit a model and simulate results 
## Apply Monte-Carlo test on fitted obj with class "list"
##############################################################
## Example 1 
## Threshold Autoregressive (TAR) Model example from TSA package
## It takes about 64.27 seconds on personal PC with 4 CPUs
library("TSA")
FitModel <- function(data){
    fit <- TSA::tar(y=log(data),p1=4,p2=4,d=3,a=0.1,b=0.9,print=FALSE)
    res <- ts(fit$std.res)
    parSpec <- list(res=res,fit=fit)
  parSpec
}
SimModel <- function(parSpec){
    fit <- parSpec$fit   
  exp(tar.sim(fit)$y)
}
data(prey.eq)
portest(FitModel(prey.eq),nslaves=4,func=list(SimModel,FitModel),pkg="TSA")
#
## Example 2
## It takes about 10.75 seconds on personal PC with 4 CPUs
FitModel <- function(data){
    fit <- ar(data,aic = FALSE, order.max=2)
    order <- 2
    res <- ts(fit$resid[-(1:order)]) 
    phi <- fit$ar
    theta <- NULL
    sigma <- fit$var.pred
    demean <- fit$x.mean
 list(res=res,phi=phi,theta=theta,sigma=sigma,demean=demean)
}
SimModel <- function(parSpec){
    res <- parSpec$res
    n <- length(res)  
    innov <- sample(x=res,size=n,replace = TRUE)  
    phi <- parSpec$phi
    theta <- parSpec$theta
    sigma <- parSpec$sigma 
    demean <- parSpec$demean 
   arima.sim(n = n, list(ar = phi, ma = theta), innov = innov, 
             sd = sqrt(sigma), mean = demean)
}
Fit <- FitModel(Nile)
portest(Fit,nslaves=4,func=list(SimModel=SimModel,FitModel=FitModel),pkg="stats")
#########################################################################
####                                                                 ####
####                Simulation using varima.sim Function             ####  
####                                                                 #### 
#########################################################################
# Simulate MA(1) where innovation series is provided via argument innov 
########################################################################
set.seed(1234)
n <- 200
phi <-  NULL
theta <-  0.6
d <- NA
sigma <- 1.9
Z <- varima.sim(phi, theta, d, sigma, n,innov=rnorm(n))
plot(Z)
########################################################################
# Simulate ARIMA(2,1,0) process with phi=c(1.3,-0.35), Gaussian innovations
# The series is truncated at lag 50
########################################################################
set.seed(1234)
Trunc.Series <- 40
n <- 1000
phi <- c(1.3, -0.35)
theta <- NULL
d <- 1
sigma <- 1
Z <- varima.sim(phi,theta,d,sigma,n,Trunc.Series=Trunc.Series)
coef(arima(Z,order=c(2,1,0)))
########################################################################
# Simulate MA(1) process with theta = 0.5, t5-distribution innovations
########################################################################
set.seed(1234)
n <- 200
phi <- NULL
theta <- 0.5
Z <- varima.sim(phi, theta, sigma=1, n=n, innov.dist="t", df=5)
plot(Z)
########################################################################
# Simulate univariate ARMA(2,1) process with length 500, 
# phi = c(1.3, -0.35), theta = 0.1. Drift equation is 8 + 0.05*t
# Stable innovations with: ALPHA = 1.75, BETA = 0, GAMMA = 1, DELTA = 0
########################################################################
set.seed(1234)
n <- 500
phi <-  c(1.3, -0.35)
theta <-  0.1
constant <- 8
trend <- 0.05
demean <- 0
d <- 0
sigma <-  0.7
ALPHA <- 1.75
BETA <- 0
GAMMA <- 1
DELTA <- 0
Stable <- c(ALPHA,BETA,GAMMA,DELTA)
Z <- varima.sim(phi,theta,d,sigma,n,constant,trend,demean,
     innov.dist="stable",StableParameters=Stable)
plot(Z)
########################################################################
# Simulate a bivariate white noise series from a multivariate t4-distribution 
########################################################################
set.seed(1234)
Z <- varima.sim(sigma=diag(2),n=200,innov.dist="t",df=4)
plot(Z)
########################################################################
# Simulate a trivariate VARMA(1,1) process with length 300. 
# phi = array(c(0.5,0.4,0.1,0.5,0,0.3,0,0,0.1), dim=c(k,k,1)), where k =3
# theta = array(c(0,0.25,0,0.5,0.1,0.4,0,0.25,0.6), dim=c(k,k,1)).
# innovations are generated from multivariate normal distribution
# The process have mean c(10, 0, 12),
# Drift equation a + b * t, where a = c(2,1,5), and b = c(0.01,0.06,0)
# The series is truncated at default value: Trunc.Series=ceiling(100/3)=34 
########################################################################
set.seed(1234)
k <- 3
n <- 300
Trunc.Series <-  50   
phi <-  array(c(0.5,0.4,0.1,0.5,0,0.3,0,0,0.1),dim=c(k,k,1))
theta <-  array(c(0,0.25,0,0.5,0.1,0.4,0,0.25,0.6),dim=c(k,k,1))
sigma <- diag(k)
constant <- c(2,1,5)
trend <- c(0.01,0.06,0)
demean <- c(10,0,12)
Z <- varima.sim(phi, theta, d = 0,sigma, n, constant,trend,demean)
plot(Z)
########################################################################
# Simulate a bivariate VARIMA(2,d,1) process with n=300, where d=(1,2). 
# phi = array(c(0.5,0.4,0.1,0.5,0,0.3,0,0),dim=c(k,k,2)),
# theta = array(c(0,0.25,0,0), dim=c(k,k,1)).
# innovations are generated from multivariate normal
# The process have mean zero and no deterministic terms.
# The variance covariance is sigma = matrix(c(1,0.71,0.71,2),2,2).
# The series is truncated at default value: Trunc.Series=ceiling(100/3)=34 
########################################################################
set.seed(1234)
k <- 2
n <- 300
Trunc.Series <-  50   
phi <-  array(c(0.5,0.4,0.1,0.5,0,0.3,0,0),dim=c(k,k,2))
theta <-  array(c(0,0.25,0,0),dim=c(k,k,1))
d <- c(1,2)
sigma <- matrix(c(1,0.71,0.71,2),k,k)
Z <- varima.sim(phi, theta, d, sigma, n)
plot(Z)
########################################################################
# Simulate a bivariate VAR(1) process with length 600. 
# Stable distribution: ALPHA=(1.3,1.6), BETA=(0,0.2), GAMMA=(1,1), DELTA=(0,0.2)
# The series is truncated at default value: Trunc.Series=min(100,200)=100 
########################################################################
set.seed(1234)
k <- 2
n <- 600
phi <- array(c(-0.2,-0.6,0.3,1.1),dim=c(k,k,1))
theta <- NULL
d <- NA
sigma <- matrix(c(1,0.71,0.71,2),k,k)
ALPHA <- c(1.3,1.6)
BETA <- c(0,0.2)
GAMMA <-c(1,1)
DELTA <-c(0,0.2)
Stable <- c(ALPHA,BETA,GAMMA,DELTA)
Z <- varima.sim(phi,theta,d,sigma,n,innov.dist="stable",StableParameters=Stable)
plot(Z)
}
}

\author{
Author: Esam Mahdi and A. Ian McLeod.

Maintainer: A. Ian McLeod <aimcleod@uwo.ca> 
}

\references{
Hosking, J. R. M. (1980). "The Multivariate Portmanteau Statistic". 
Journal of American Statistical Association, 75, 602-608.

Li, W. K. and McLeod, A. I. (1981). "Distribution of The Residual Autocorrelations in 
Multivariate ARMA Time Series Models". 
Journal of The Royal Statistical Society, Series B, 43, 231-239.

Lin, J.-W. and McLeod, A.I. (2006). "Improved Generalized Variance Portmanteau Test".
Computational Statistics and Data Analysis 51, 1731-1738.

Lin, J.-W. and McLeod, A.I. (2008). "Portmanteau Tests for ARMA Models with Infinite Variance".  
Journal of Time Series Analysis, 29, 600-617.

Mahdi, E. and McLeod, A.I. (2012). "Improved Multivariate Portmanteau Test". 
Journal of Time Series Analysis, 33(2), 211-222.

McCulloch, J. H. (1986). "Simple Consistent Estimator of Stable Distribution Parameters".
Commun. Statist.--Simula., 15(4), 1109-1136.

McLeod A.I, Li W.K (1983). "Distribution of the Residual Autocorrelation in Multivariate
ARMA Time Series Models". Journal of Time Series Analysis, 4, 269-273.

McLeod, A.I., Yu, Hao, and Krougly, Z. L. (2007). "Algorithms for Linear Time Series Analysis".
Journal of Statistical Software.

Pena, D. and Rodriguez, J. (2006). "The log of the determinant of the autocorrelation
matrix for testing goodness of fit in time series". Journal of Statistical Planning and
Inference, 136, 2706-2718.

Tierney, L., Rossini, A. J., Li, N., and Sevcikova, H. (2009). snow: Simple Network of Workstations.
\code{R} package version 0.3-10. \url{http://CRAN.R-project.org/package=snow}.

Wuertz, D. and core team members R (2012). fGarch: Rmetrics - Autoregressive Conditional
Heteroskedastic Modelling. \code{R} package version 2150.81. \url{http://CRAN.R-project.org/package=fGarch}.

Yu, H. (2002). Rmpi: Parallel Statistical Computing in R. R News, 2(2), 10-14. 
\url{http://CRAN.R-project.org/doc/Rnews}.


}

\keyword{ PACKAGE }
\keyword{ ts}




