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
     NREP=1000,InfiniteVarianceQ=FALSE,SquaredQ=FALSE,SetSeed=TRUE)
}

\arguments{
  \item{obj}{if \code{obj} is an object of class \code{"ar"}, \code{"arima0"}, \code{"Arima"}, \code{"varest"}, 
             \code{"FitAR"}, or \code{"FitFGN"} then a portmanteau goodness-of-fit test is done on the fitted model. 
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
  \item{SquaredQ}{if \code{TRUE} then apply the test on the squared values. 
             This checks for Autoregressive Conditional Heteroscedastic, \code{ARCH}, effects. 
             When \code{SquaredQ = FALSE}, then apply the test on the usual residuals.} 
  \item{SetSeed}{if \code{TRUE} then \code{set.seed} is initialized.} 
}

\details{
The portmanteau test statistics, \code{\link{gvtest}}, \code{\link{BoxPierce}}, \code{\link{LjungBox}},
\code{\link{Hosking}}, and \code{\link{LiMcLeod}} are implemented based on the Monte-Carlo techniques
and the approximation asymptotic distributions. 
The null hypothesis assuming that the fitted model is an adequate 
model and the residuals behave like white noise series. 
The \code{snow} package must be installed if one 
decide to choose the argument \code{MonteCarlo=TRUE} provided that 
\code{nslaves>1}.
This function can be used for testing the  
adequacy in the nonseasonal fitted \code{ARIMA}, \code{VAR}, 
and Fractional Gaussian Noise, \code{FGN}, models.
Also, it can be used to check for randomness as well as for \code{ARCH} effects.
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

Mahdi, E. and McLeod, A.I. (2011). "Improved Multivariate Portmanteau Diagnostic Test". Submitted.

McLeod A.I, Li W.K (1983). "Distribution of the Residual Autocorrelation in Multivariate
ARMA Time Series Models". Journal of Time Series Analysis, 4, 269-273.

Pena, D. and Rodriguez, J. (2006). "The log of the determinant of the autocorrelation
matrix for testing goodness of fit in time series". Journal of Statistical Planning and
Inference, 136, 2706-2718.
}

\seealso{
 \code{\link{varima.sim}}, \code{\link{ar}}, \code{\link{arima0}}, \code{\link{arima}}, 
\code{\link[forecast]{Arima}}, \code{\link[forecast]{auto.arima}}, 
\code{\link{FitAR}}, \code{vars}, \code{FitFGN}, \code{\link{BoxPierce}},
\code{\link{gvtest}}, \code{\link{LjungBox}},
 \code{\link{Hosking}}, \code{\link{LiMcLeod}}, \code{\link{fitstable}}
 }

\examples{
###############################################################
## Example 1
###############################################################
# Simple example using Monte-Carlotest in portes.
# 100 replications takes about 6 seconds on single cpu 2.67 GHz.
##
ans <- arima(Nile, order=c(1,0,1))
portest(ans, NREP=100)
##
###############################################################
## Example 2
###############################################################
# Simulate a bivariate VARIMA(1,c(0,0),1) process with length 300. 
# Apply gvtest based on the two methods implemented in portes. 
# 100 replications takes about 23.63 seconds on single cpu 2.67 GHz.
##
k <- 2
n <- 300
Trunc.Series <-  50   
phi <-  array(c(0.5,0.4,0.1,0.5),dim=c(k,k,1))
theta <-  array(c(0,0.25,0,0),dim=c(k,k,1))
d <- NA
sigma <- matrix(c(1,0.71,0.71,2),k,k)
z <- varima.sim(phi,theta,d,sigma,n)
ans <- ar(z)                  
portest(ans, MonteCarlo=FALSE)  ## asymptotic distribution method
portest(ans, NREP=100)          ## Monte-Carlo method
##
###############################################################
## Example 3
###############################################################
# Simulate a bivariate VARIMA(2,d,0) process with length 300.
## d is a vector (1,2)
## Fit VAR(2) using the function VAR in the package "vars" 
# Apply gvtest based on the two methods implemented in portes. 
# 100 replications takes about 30.17 seconds on single cpu 2.67 GHz.
##
\dontrun{
library(vars)
k <- 2
n <- 300
Trunc.Series <-  50   
phi <-  array(c(0.5,0.4,0.1,0.5,0,0.25,0,0),dim=c(k,k,2))
theta <-  NULL
d <- c(1,2)
sigma <- matrix(c(1,0.71,0.71,2),k,k)
z <- varima.sim(phi,theta,d,sigma,n)
ans <- VAR(z, p=2)            ## inadequate fitted model
portest(ans, MonteCarlo=FALSE) ## asymptotic distribution method
portest(ans, NREP=100)         ## Monte-Carlo method
##
###############################################################
## Example 4
###############################################################
## Checks the residuals for randomness using LjungBox test. 
portest(rnorm(100),test="LjungBox")
##
###############################################################
## Example 5 - Using "snow" package 
###############################################################
## Apply Hosking test on fitted VAR(2) model to WestGerman data. 
##
library("snow")
data("WestGerman")
DiffData <- matrix(numeric(3 * 91), ncol = 3)
for (i in 1:3) DiffData[, i] <- diff(log(WestGerman[, i]), lag = 1)
Fit2 <- ar.ols(DiffData, aic=FALSE, order.max = 2, intercept = FALSE)
portest(Fit2,test="Hosking",nslaves=8)        ## Monte-Carlo
portest(Fit2,test="Hosking",MonteCarlo=FALSE) ## asymptotic distribution                       
##
###############################################################
## Example 6 - Using "snow" package 
###############################################################
## Test monthly log stock returns of Intel data for ARCH effects.
## gvtest statistic on PC with 8 CPU's using "snow". 
## It takes 16.75 seconds based on the Monte-Carlo test.
##
data(monthintel)
returns <- as.ts(monthintel)
lags <- c(10, 20, 30, 40)
portest(returns,lags=lags,MonteCarlo=TRUE,nslaves=8,SquaredQ=TRUE) 
##
###############################################################
## Example 7 - Using "snow" package 
###############################################################
## Fit Fractional Gaussian Noise, FGN, to NileMin data in FGN package.
## Monte-Carlo of gvtest on 8 CPU's using "snow".
##
library(FGN)
data(NileMin)
NILE.FGN <- FitFGN(NileMin)
lags <- c(5, 10, 20)
##				
## gvtest statistic on fitted model (55 seconds)
portest(NILE.FGN, lags=lags, nslaves=8) 
##
## gvtest statistic on residuals (6 seconds)
res <- NILE.FGN$res
portest(res, lags=lags, nslaves=8)
##
}
}

\keyword{htest}
\keyword{ ts}


