\name{portest}
\alias{portest}

\title{Portmanteau Test Statistics}

\description{
 Univariate or multivariate portmanteau 
test statistics based on the Monte-Carlo techniques or asymptotic distributions.
}

\usage{
  portest(obj,lags=seq(5,30,5),order=0,test=c("gvtest","BoxPierce",
    "LjungBox","Hosking","LiMcLeod"),MonteCarlo=TRUE,nslaves=1,
     NREP=1000,InfiniteVarianceQ=FALSE,SquaredQ=FALSE,
     func=list(SimModel=NULL,FitModel=NULL),pkg=NULL,SetSeed=TRUE)
}

\arguments{
  \item{obj}{if \code{obj} is an object of class \code{"ar"}, \code{"arima0"}, \code{"Arima"}, \code{"varest"}, 
             \code{"FitAR"}, \code{"FitFGN"}, \code{"garch"}, \code{"fGARCH"}, or \code{"list"} then a portmanteau 
             goodness-of-fit test is done on the fitted model. 
             Otherwise, for \code{obj} with class \code{"ts"}, \code{"numeric"}, \code{"matrix"}, 
             or \code{("mts" "ts")}, a test of randomness is done.}
  \item{lags}{vector of lag values is used for portmanteau test.}
  \item{order}{as described in \code{BoxPierce}, \code{gvtest}, \code{Hosking}, \code{LiMcLeod}, and \code{LjungBox}
             and needed only when \code{MonteCarlo = FALSE} is selected. 
             }   
  \item{test}{portmanteau test statistic type}
  \item{MonteCarlo}{if \code{TRUE} then apply the Monte-Carlo version of portmanteau statistic. 
             Otherwise, apply the asymptotic distribution. }
  \item{nslaves}{number of slaves needed to use in parallel calculations.
             Default is one single CPU.}
  \item{NREP}{number of replications needed for Monte-Carlo test.}
  \item{InfiniteVarianceQ}{\code{FALSE}, assumes innovations with finite variance. 
             Otherwise, innovations follow stable distribution with infinite variance.}
  \item{SquaredQ}{as described in \code{BoxPierce}, \code{gvtest}, 
             \code{Hosking}, \code{LiMcLeod}, and \code{LjungBox}.} 
  \item{func}{additional argument defined as a list with two specified functions, 
             \code{SimModel} and \code{FitModel}.
             This argument is needed when the class of \code{obj} is \code{"list"} 
             (see details and following example).} 
  \item{pkg}{the name of the library to be loaded if the Monte-Carlo significance test is used with
            an object \code{obj} with class \code{"list"}.} 
  \item{SetSeed}{if \code{TRUE} then \code{set.seed} is initialized.} 
}

\details{
The portmanteau test statistics, \code{\link{gvtest}}, \code{\link{BoxPierce}}, \code{\link{LjungBox}},
\code{\link{Hosking}}, and \code{\link{LiMcLeod}} are implemented based on the Monte-Carlo techniques
and the approximation asymptotic distributions as described in Mahdi and McLeod (2012). 
The null hypothesis assuming that the fitted model is an adequate 
model and the residuals behave like white noise series. 
The parallel computing program \code{MPICH2} must be properly installed
in the I O environment system in which the \code{portes} package will run if one 
decide to choose the argument \code{MonteCarlo=TRUE} provided that 
\code{nslaves} equals to a positive integer more than 1.
This function can be used for testing the  
adequacy in the nonseasonal fitted \code{ARIMA}, \code{VAR}, 
and Fractional Gaussian Noise, \code{FGN}, models.
Also, it can be used to check for randomness as well as for \code{GARCH} effects.
Any other fitted model, for example, threashold autoregression model,
may also be tested for adequacy.
In this case, two functions, \code{SimModel()} and \code{FitModel()},
must be provided via the argument \code{func}.
The object \code{obj} is the output of the fitted model coded in
the function \code{FitModel} and it is a \code{"list"} with at least
\code{res}, the residuals from the fitted model in \code{FitModel()}.
The output from the function \code{SimModel()} is a simulated 
univariate or multivariate series from the fitted model coded in the function
\code{FitModel()}.
The argument \code{pkg} represents the name of the \code{R} package where the fitted
model build in (see the last two examples). 
}

\value{
The portmanteau test statistic with the associated p-values 
for different lag values.
}

\author{ Esam Mahdi and A.I. McLeod.}

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

\seealso{
 \code{\link{varima.sim}}, \code{\link{ar}}, \code{\link{arima0}}, \code{\link{arima}}, 
\code{\link[forecast]{Arima}}, \code{\link[forecast]{auto.arima}}, 
\code{\link[FitAR]{FitAR}}, \code{\link[vars]{VAR}}, \code{\link[FGN]{FitFGN}}, \code{\link{BoxPierce}},
\code{\link{gvtest}}, \code{\link{LjungBox}},
 \code{\link{Hosking}}, \code{\link{LiMcLeod}}, \code{\link{fitstable}},
\code{\link[tseries]{garch}}, \code{\link[fGarch]{garchFit}}, 
\code{\link[TSA]{tar}}.
 }

\examples{
#to save time the other examples are not run
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
}
}

\keyword{htest}



