\name{portes-package}
\alias{portes-package}
\docType{package}

\title{
Portmanteau Tests for Univariate and Multivariate Time Series Models 
}

\description{
This package contains a set of portmanteau diagnostic checks for univariate and multivariate time series
based on the asymptotic approximation distributions and the Monte-Carlo significance test.
More details about the portmanteau test statistics are available online from the vignette of this package.
This package can be also used to simulate univariate and multivariate data 
from seasonal/nonseasonal \code{ARIMA}/ \code{VARIMA} models with
innovations from finite or infinite variances from stable distributions. 
The simulated data may have deterministic terms, constant drifts and time trends, with non-zero means.
}

\details{
\tabular{ll}{
Package: \tab portes\cr
Type: \tab Package\cr
Version: \tab 2.1-4\cr
Date: \tab 2017-05-05\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
Depends: \tab \code{R (>= 3.2.5), parallel(V.3.2.1)}\cr 
Suggests: \tab \code{akima(V.0.6-2), fGarch(V.3010.82.1), FGN(V.2.0-12),}\cr
\tab FitAR(V.1.94), \code{forecast(V.8.0), fracdiff(V.1.4-2), gstat(V.1.1-5),}\cr
\tab TSA(V.1.01), tseries(V.0.10-38), \code{vars(V.1.5-2)}\cr
Classification/ACM: \tab G.3, G.4, I.5.1\cr
Classification/MSC: \tab 62M10, 91B84\cr
License: \tab GPL (>= 2)\cr
}
}

\section{Main Function}{

The main function in this package, \code{\link{portest}}, 
is used with univariate and multivariate time series.
It implements any possible test statistic including the Portmanteau 
test statistics, \code{\link{MahdiMcLeod}}, 
\code{\link{BoxPierce}}, \code{\link{LjungBox}}, \code{\link{Hosking}}, 
\code{\link{LiMcLeod}}, and others based on two methods.
The first method uses the Monte-Carlo techniques as described in 
Lin and McLeod (2006) and Mahdi and McLeod (2012), whereas the second one uses 
the approximation asymptotic distribution.

Originally, the generalized variance portmanteau test of \code{\link{MahdiMcLeod}} 
for univariate time series was derived by Pena and Rodriguez (2002)
based on the gamma distribution. 
Lin and McLeod (2006) proposed the Monte-Carlo version of this test and 
Mahdi and McLeod (2012) extended both methods to the multivariate case. 
Simulation results suggest that the Monte-Carlo version of 
\code{\link{MahdiMcLeod}} statistic is more accurate
and powerful than its competitors proposed by Box and Pierce (1970), Ljung and Box (1978),
and Pena and Rodriguez (2002, 2006) in the univariate time series as well as 
Hosking (1980) and Li and McLeod (1981) in the multivariate time series.

The powerful parallel computing framework facility, using the package \code{parallel} 
(Gonzalo Vera, Ritsert Jansen and Remo Suppi 2008),
is implemented in this function.
The \code{parallel} package handles running much larger chunks of computations in parallel and was first 
included in \code{R 2.14.0} based on the work done for \code{CRAN} packages
\code{multicore} and \code{snow}.

The default argument in \code{\link{portest}} function, \code{ncores=1}, 
implements the Monte-Carlo test on PC with only one CPU. 
Set the argument \code{ncores} equals to a positive integer number greater than 1, 
provided that the default argument \code{MonteCarlo=TRUE} is selected, 
the package \code{parallel} will be implemented for faster calculations.
When \code{MonteCarlo=FALSE} is selected, the test statistic 
selected from the argument \code{test} will be implemented
based on the asymptotic approximation distribution.
The default test statistic is the generalized variance test, \code{\link{MahdiMcLeod}}.

\subsection{Test for usual residuals and \code{ARCH} effects}{

By setting the argument \code{squared.residuals=TRUE} in \code{\link{portest}} function,
the portmanteau test using the asymptotic distribution approximation or the 
Monte-Carlo significance test 
(depending on the choice of the argument \code{MonteCarlo} whether \code{FALSE} or \code{TRUE})
for \code{ARCH} effects will be implemented on the squared residuals.
Otherwise, the portmanteau test will be applied on the usual residuals when the 
default argument \code{squared.residuals=FALSE} is selected.
}

\subsection{Test for seasonality}{ 
This package is also useful in testing for seasonality
as discussed by McLeod (1978). 
The Monte-Carlo significance tests available from the function \code{\link{portest}} can be used to test for 
the seasonality in the univariate \code{SARIMA} models fitted by the functions \code{Arima()} 
or \code{auto.arima()} available from the package \code{forecast}. 
With the asymptotic distribution of the portmanteau tests, 
the seasonal period is entered the argument \code{season},
where \code{season = 1} is used for usual test with no seasonality check.
}

\subsection{Monte-Carlo test for residuals with infinite variances (heavy tails in stable distributions)}{

By selecting the argument \code{innov.dist = stable}, the portmanteau test statistics on residuals with
infinite variances from stable processes are implemented.
}

\subsection{Test for adequacy of fitted time series regression models}{

The \code{\link{portest}} function has the ability to test the adequacy of fitted 
linear models, \code{LM} and generalized linear models, \code{GLM} that has time series disturbances.
It computes the test statistics \code{\link{MahdiMcLeod}}, \code{\link{BoxPierce}}, \code{\link{LjungBox}}, 
\code{\link{Hosking}}, and \code{\link{LiMcLeod}} portmanteau tests,
based on both methods: the asymptotic distribution and the Monte-Carlo significance test.   
This function also involves testing the adequacy of the fitted
linear models using the Monte-Carlo version of the univariate \code{Generalized Durbin-Watson} test statistic.
This can be done by selecting the arguments \code{test = "other"} and \code{fn = dwt}, where \code{dwt}
is a coded written \code{R} function calculates the \code{Generalized Durbin-Watson} test statistic.
}

\subsection{Goodness-of-fit test for any fitted model}{

The portmanteau test statistics implemented in the functions, \code{\link{MahdiMcLeod}}, \code{\link{BoxPierce}}, 
\code{\link{LjungBox}}, \code{\link{Hosking}}, \code{\link{LiMcLeod}}, can be 
used for testing the adequacy of any fitted model, such as \code{ARFIMA}, \code{ARCH}, \code{GARCH}, 
and \code{TAR} models, using the Monte-Carlo significance test 
as well as the asymptotic approximation distribution test. 
Note that the argument function \code{fn} may used with the argument \code{test = "other"} 
to implement any test statistic other than those mentioned in the previous list.
More details with plenty illustrative examples are given in the documentation of
the main function \code{\link{portest}}.
}
}

\section{Simulate data from seasonal/nonseasonal \code{ARIMA}/ \code{VARIMA}}{

The function \code{\link{varima.sim}} in this package is useful for simulating data
from seasonal/nonseasonal \code{ARIMA}/ \code{VARIMA} of order \eqn{(p,d,q)\times(ps,ds,qs)_s} 
with or without deterministic terms (drift and trend).
The innovations series can be inserted in the function \code{portest()} 
via its argument \code{innov} as an initial series to generate data using nonparametric 
bootstrap procedure or stable distribution specified from the argument \code{innov.dist}.
The arguments \code{d} and \code{D} must entered as a nonnegative integer 
in the \code{ARIMA} and \code{SARIMA} cases, whereas it must entered as a vector of \eqn{k} components 
\eqn{d_1,...,d_k} and \eqn{ds_1,...,ds_k} in the \code{VARIMA} and \code{seasonal SVARIMA} cases.
\eqn{d_i} represents the usual difference lag need to be applied on series \eqn{i} and
\eqn{ds_i} represents the seasonal difference lag need to be applied on series \eqn{i}.
The components of the argument \code{par.stable} are the 
stable parameters \code{Alpha, Beta, Scale}, and \code{Location} are needed to generate  
innovations from stable distribution. 
}

\examples{
\dontrun{
#################################################################################
####                                                                         ####
####                             Portmanteau Tests                           ####  
####                                                                         #### 
#################################################################################
## Monte-Carlo (MC) and asymptotic tests for randomness series                 ##
#################################################################################
data("DEXCAUS")
returns <- log(DEXCAUS[-1]/DEXCAUS[-length(DEXCAUS)])
portest(returns)                    ## MC using one CPU takes about 24.16 seconds
portest(returns, ncores=4)          ## MC using 4 CPUs takes about 9.51 seconds
portest(returns, MonteCarlo=FALSE)                      ## asymptotic MahdiMcLeod
portest(returns,test="LjungBox", MonteCarlo=FALSE)      ## asymptotic LjungBox 
#################################################################################
## Monte-Carlo goodness-of-fit arima test using 4 CPUs                         ##
#################################################################################
## arima() function takes about 11.32 seconds
## Example 1
ans1 <- arima(WWWusage,order=c(3,1,0))
portest(ans1, ncores = 4)
#
## arima0() function takes about 11.1 seconds
## Example 2
ans2 <- arima0(WWWusage,order=c(3,1,0))
portest(ans2, ncores = 4)
#
## Arima() or auto.Arima() function takes about 12.1 seconds
## Example 3
require("forecast")
ans3 <- Arima(WWWusage,order=c(3,1,0))
portest(ans3, ncores = 4)
detach(package:forecast)
#
## ar() function takes about 7.39 seconds
## Example 4
ans4 <- ar(Nile,order.max=2)
portest(ans4, ncores = 4)
#
## ar() function with your own R code takes about 8.75 seconds
## Example 5
fit.model <- function(data){
    fit <- ar(data,aic = FALSE, order.max=2)
    order <- 2
    res <- ts(fit$resid[-(1:order)]) 
    phi <- fit$ar
    theta <- NULL
    sigma <- fit$var.pred
    demean <- fit$x.mean
 list(res=res,phi=phi,theta=theta,sigma=sigma,demean=demean)
}
sim.model <- function(parSpec){
    res <- parSpec$res
    n <- length(res)  
    innov <- function(n) ts(stats::rnorm(n, mean = demean, sd = sqrt(sigma))) 
    phi <- parSpec$phi
    theta <- parSpec$theta
    sigma <- parSpec$sigma 
    demean <- parSpec$demean 
   arima.sim(n = n, list(ar = phi, ma = theta), rand.gen=innov)
}
ans5 <- fit.model(Nile)
portest(ans5,ncores=4,model=list(sim.model=sim.model,fit.model=fit.model),pkg.name="stats")
#################################################################################
## Monte-Carlo test for seasonality                                            ##
#################################################################################
## Accidental Deaths in the US 1973 - 1978
require("forecast")
seasonal.arima<-Arima(USAccDeaths,order=c(0,1,1),seasonal=list(order= c(0,1,1)))
portest(seasonal.arima,ncores=4,nrep=1000,lags=1:5)
## Quarterly U.K. economic time series from 1957 Q3 to 1967 Q4
cd <- EconomicUK[,1]
cd.fit <- Arima(cd,order=c(0,1,0),seasonal=list(order=c(0,1,1),period=4))
portest(cd.fit, lags = c(5,10),ncores=4)
detach(package:forecast)
#################################################################################
## Monte-Carlo test for linear models and time series regression               ##
#################################################################################
## Linear Model
require("car")
fit <- lm(fconvict ~ tfr + partic + degrees + mconvict, data=Hartnagel)
portest(fit,lags=1:3,ncores=4)          ## MC of MahdiMcLeod test
## MC of generalized Durbin-Watson test needs the argument function fn() as follows
fn <- function(obj,lags){
     test.stat <- numeric(length(lags))
       for (i in 1:length(lags))
          test.stat[i] <- -sum(diff(obj,lag=lags[i])^2)/sum(obj^2)
       test.stat
}
portest(fit,lags=1:3,test="other",fn=fn,ncores=4) 
detach(package:car)
## Time series regression 
require("forecast")
fit.arima <- Arima(LakeHuron, order = c(2,0,0), xreg = time(LakeHuron)-1920)
portest(fit.arima,ncores=4)
detach(package:forecast)
#################################################################################
## Monte-Carlo goodness-of-fit VAR test - Multivariate series                  ##
#################################################################################
data("IbmSp500")
ibm <- log(IbmSp500[,2]+1)*100
sp500 <- log(IbmSp500[,3]+1)*100
IBMSP500 <- data.frame(cbind(ibm,sp500))
## ar.ols() function takes about 9.11 seconds
ans6 <-  ar.ols(IBMSP500, aic=FALSE, intercept=TRUE, order.max=5)
portest(ans6,nrep=100,test="MahdiMcLeod",ncores=4,innov.dist="t",dft=5)
## VAR() function takes about 11.55 seconds
require("vars")
ans7 <- VAR(IBMSP500, p=5) 
portest(ans7,nrep=100,test="MahdiMcLeod",ncores=4,innov.dist="bootstrap")
portest(ans7,test="Hosking",MonteCarlo=FALSE)          ## asymptotic Hosking test
detach(package:vars)
#################################################################################
## Monte-Carlo test for models with heavy tails stable distributions           ##
#################################################################################
## It takes about 32.7 seconds on personal PC with 4 CPUs 
data("CRSP")
CRSP.AR5<- arima(CRSP, c(5, 0, 0))
lags <- c(10, 20, 30)
portest(CRSP.AR5,lags=lags,ncores=4,nrep=1000,innov.dist="stable")   
#################################################################################
## Monte-Carlo test for ARCH/GARCH effects using 4 CPUs                        ##
#################################################################################
## It takes about 12.65 seconds
data("monthintel")
returns <- as.ts(monthintel)
lags <- c(5, 10, 20, 40)
portest(returns, lags = lags, ncores = 4, squared.residuals = FALSE)
portest(returns,lags=lags,ncores=4,squared.residuals=TRUE,innov.dist="t",dft=5)
#################################################################################
## Monte-Carlo test for Threshold Autoregressive (TAR) Models                  ##
## It takes about 54.27 seconds on personal PC with 4 CPUs                     ##
#################################################################################
require("TSA")
fit.model <- function(data){
    fit <- TSA::tar(y=log(data),p1=4,p2=4,d=3,a=0.1,b=0.9,print=FALSE)
    res <- ts(fit$std.res)
    parSpec <- list(res=res,fit=fit)
  parSpec
}
sim.model <- function(parSpec){
    fit <- parSpec$fit   
  exp(tar.sim(fit)$y)
}
data(prey.eq)
portest(fit.model(prey.eq),ncores=4,model=list(sim.model,fit.model),pkg.name="TSA")
detach(package:TSA)
#################################################################################
####                                                                         ####
####                   Simulation using varima.sim Function                  ####  
####                                                                         #### 
#################################################################################
# Simulate white noise series from a Gaussian distribution                      #
#################################################################################
set.seed(1234)
Z1 <- varima.sim(n=400)       ## a univariate series
plot(Z1)
Z2 <- varima.sim(n=400,k=2)   ## a bivariate series
plot(Z2)
Z3 <- varima.sim(n=400,k=5)   ## a multivariate series of dimension 5
plot(Z3)
#################################################################################
# Simulate MA(1) where innovation series is provided via argument innov         #
#################################################################################
set.seed(1234)
n <- 200
theta <-  0.6
Z<-varima.sim(list(ma=theta),n=n,innov=rnorm(n),innov.dist="bootstrap")
plot(Z)
#################################################################################
# Simulate seasonal ARIMA(2,1,0)*(0,2,1)_12 process with phi=c(1.3,-0.35),      #
# theta.season = 0.8. Gaussian innovations. The series is truncated at lag 50   #
#################################################################################
set.seed(12834)
n <- 100
phi <- c(1.3, -0.35)
theta.season <- 0.8
Z<-varima.sim(list(ar=phi,d=1,sma=theta.season,D=2),n=n,trunc.lag=50)
plot(Z)
acf(Z)
#################################################################################
# Simulate seasonal ARMA(0,0,0)*(2,0,0)_4 process with intercept = 0.8          #
# phi.season = c(0.9,-0.9), period = 4, t5-distribution innovations with df = 3 #
#################################################################################
set.seed(1234)
n <- 200
phi.season <- c(0.9,-0.9)
Z<-varima.sim(list(sar=phi.season,period=4),n=n,constant=0.8,innov.dist="t",dft=3)
plot(Z)
acf(Z)
arima(Z,order=c(0,0,0),seasonal = list(order = c(2,0,0),period=4))
#################################################################################
# Simulate univariate ARMA(2,1) process with length 500,                        #
# phi = c(1.3, -0.35), theta = 0.1. Drift equation is 8 + 0.05*t                #
# Stable innovations with: Alpha = 1.75, Beta = 0, Scale = 1, Location = 0      #
#################################################################################
set.seed(1234)
n <- 500
phi <-  c(1.3, -0.35)
theta <-  0.1
Alpha <- 1.75
Beta <- 0
Scale <- 1
Location <- 0
Stable <- c(Alpha,Beta,Scale,Location)
Z <- varima.sim(list(ar=phi,ma=theta),n=n,constant=8,trend=0.05,demean=0,
     innov.dist="stable",par.stable=Stable)
plot(Z)
#################################################################################
# Simulate a bivariate white noise series from a multivariate t4-distribution   #
# Then use the nonparametric bootstrap method to generate a seasonal SVARIMA    #
# of order (0,d,0)*(0,0,1)_12 with d = c(1, 0), n= 250, k = 2, and              #
# theta.season=array(c(0.5,0.4,0.1,0.3),dim=c(k,k,1))                           #
#################################################################################
set.seed(1234)
Z1 <- varima.sim(n=250,k=2,innov.dist="t",dft=4)
theta.season=array(c(0.5,0.4,0.1,0.3),dim=c(2,2,1)) 
Z2 <- varima.sim(list(sma=theta.season,d=c(1,0)),n=250,k=2,
                 innov=Z1,innov.dist="bootstrap")
plot(Z2)
#################################################################################
# Simulate a bivariate VARIMA(2,d,1) process with length 300, where d=(1,2).    #
# phi = array(c(0.5,0.4,0.1,0.5,0,0.3,0,0),dim=c(k,k,2)),                       #
# theta = array(c(0,0.25,0,0), dim=c(k,k,1)).                                   #
# innovations are generated from multivariate normal                            #
# The process have mean zero and no deterministic terms.                        #
# The variance covariance is sigma = matrix(c(1,0.71,0.71,2),2,2).              #
# The series is truncated at default value: trunc.lag=ceiling(100/3)=34         #
#################################################################################
set.seed(1234)
k <- 2
n <- 300 
phi <-  array(c(0.5,0.4,0.1,0.5,0,0.3,0,0),dim=c(k,k,2))
theta <-  array(c(0,0.25,0,0),dim=c(k,k,1))
d <- c(1,2)
sigma <- matrix(c(1,0.71,0.71,2),k,k)
Z <- varima.sim(list(ma=phi,ma=theta,d=d),n=n,k=2,sigma=sigma)
plot(Z)
#################################################################################
# Simulate a trivariate Vector SVARMA(1,0,0)*(1,0,0)_12 process with length 300 #
# phi = array(c(0.5,0.4,0.1,0.5,0,0.3,0,0,0.1), dim=c(k,k,1)), where k =3       #
# phi.season = array(c(0,0.25,0,0.5,0.1,0.4,0,0.25,0.6), dim=c(k,k,1)).         #
# innovations are generated from multivariate normal distribution               #
# The process have mean c(10, 0, 12),                                           #
# Drift equation a + b * t, where a = c(2,1,5), and b = c(0.01,0.06,0)          #
# The series is truncated at default value: trunc.lag=ceiling(100/3)=34         #
#################################################################################
set.seed(1234)
k <- 3
n <- 300
phi <-  array(c(0.5,0.4,0.1,0.5,0,0.3,0,0,0.1),dim=c(k,k,1))
phi.season <-  array(c(0,0.25,0,0.5,0.1,0.4,0,0.25,0.6),dim=c(k,k,1))
constant <- c(2,1,5)
trend <- c(0.01,0.06,0)
demean <- c(10,0,12)
Z <- varima.sim(list(ar=phi,sar=phi.season),n=n,k=3,constant=constant,
trend=trend,demean=demean)
plot(Z)
acf(Z)
#################################################################################
# Simulate a bivariate VAR(1) process with length 600.                          #
# Stable dist.: Alpha=(1.3,1.6), Beta=(0,0.2), Scale=(1,1), Location=(0,0.2)    #
# The series is truncated at default value: trunc.lag=min(100,200)=100          #
#################################################################################
set.seed(1234)
k <- 2
n <- 600
phi <- array(c(-0.2,-0.6,0.3,1.1),dim=c(k,k,1))
theta <- array(c(1,-0.2,0.71,0.2),dim=c(k,k,1))
Alpha <- c(1.3,1.6)
Beta <- c(0,0.2)
Scale <-c(1,1)
Location <-c(0,0.2)
Stable <- c(Alpha,Beta,Scale,Location)
Z<-varima.sim(list(ar=phi,ma=theta),n=n,k=2,innov.dist="stable",par.stable=Stable)
plot(Z)
}
}

\author{
Author: Esam Mahdi and A. Ian McLeod.

Maintainer: A. Ian McLeod <aim@stats.uwo.ca>
}

\references{
Box, G.E.P. and Pierce, D.A. (1970). "Distribution of Residual Autocorrelation in Autoregressive-Integrated Moving Average
Time Series Models". Journal of American Statistical Association, 65, 1509-1526.

Hosking, J. R. M. (1980). "The Multivariate Portmanteau Statistic". 
Journal of American Statistical Association, 75, 602-608.

Li, W. K. and McLeod, A. I. (1981). "Distribution of The Residual Autocorrelations in 
Multivariate ARMA Time Series Models". Journal of The Royal Statistical Society, Series B, 43, 231-239.

Lin, J.-W. and McLeod, A.I. (2006). "Improved Generalized Variance Portmanteau Test".
Computational Statistics and Data Analysis 51, 1731-1738.

Lin, J.-W. and McLeod, A.I. (2008). "Portmanteau Tests for ARMA Models with Infinite Variance".  
Journal of Time Series Analysis, 29, 600-617.

Ljung, G.M. and Box, G.E.P (1978). "On a Measure of Lack of Fit in Time Series Models". 
Biometrika, 65, 297-303.

Mahdi, E. and McLeod, A.I. (2012). "Improved Multivariate Portmanteau Test". 
Journal of Time Series Analysis, 33(2), 211-222.

McLeod A.I, Li W.K (1983). "Distribution of the Residual Autocorrelation in Multivariate
ARMA Time Series Models". Journal of Time Series Analysis, 4, 269-273.

McLeod A.I, (1978). On the distribution and applications of residual autocorrelations in
Box-Jenkins modelling. Journal of the Royal Statistical Society B, 40(3), 296-302.

Pena, D. and Rodriguez, J. (2002). "A Powerful Portmanteau Test of Lack of Test for Time Series".
Journal of American Statistical Association, 97, 601-610.

Pena, D. and Rodriguez, J. (2006). "The log of the determinant of the autocorrelation
matrix for testing goodness of fit in time series". Journal of Statistical Planning and
Inference, 136, 2706-2718.

Tierney, L., Rossini, A. J., Li, N., and Sevcikova, H. (2016). snow: Simple Network of Workstations.
\code{R} package version 0.4-2. \url{https://CRAN.R-project.org/package=snow}.

Gonzalo Vera and Ritsert C. Jansen and Remo L. Suppi (2008). R/parallel - speeding 
up bioinformatics analysis with R. BMC Bioinformatics, 9:390.
}

\keyword{ PACKAGE }
\keyword{ ts}




