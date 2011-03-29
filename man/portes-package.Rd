\name{portes-package}
\alias{portes-package}
\docType{package}

\title{
Portmanteau Tests for Time Series Models 
}

\description{
This package contains a set of portmanteau diagnostic checks for univariate and multivariate time series
based on the asymptotic approximation distributions and the Monte-Carlo significance test.
It can be used for simulating data from nonseasonal \code{ARIMA} or \code{VARIMA} Models whith
innovations from finite or infinite variances distributions. 
The simulated data may have deterministic terms, a constant drift and a time trend, with non-zero mean.
}

\details{
\tabular{ll}{
Package: \tab portes\cr
Type: \tab Package\cr
Version: \tab 1.07\cr
Date: \tab 2011-03-29\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
Depends: \tab FitAR, akima \cr
Suggests: \tab FGN, forecast, vars, snow \cr
Recommends: \tab rlecuyer \cr
Classification/ACM: G.3, G.4, I.5.1\cr
Classification/MSC: 62M10, 91B84\cr
License: GPL (>= 2)\cr
}
}

\section{Main Function}{
The main function in this package, \code{\link{portest}}, 
is used with univariate and multivariate time series.
It implements the Portmanteau test statistics, \code{\link{gvtest}}, 
\code{\link{BoxPierce}}, \code{\link{LjungBox}}, \code{\link{Hosking}}, 
and \code{\link{LiMcLeod}} based on two methods.
The first method uses the Monte-Carlo techniques as described by 
Lin and McLeod (2006), Mahdi and McLeod (2011) and the second one uses 
the approximation asymptotic distribution.
Originally, the generalized variance portmanteau test, \code{\link{gvtest}}, for univariate time series 
based on the gamma distribution was proposed by Pena and Rodriguez (2002). 
Lin and McLeod (2006) proposed the Monte-Carlo version of this test and 
Mahdi and McLeod (2011) extended both methods to the multivariate case. 
Simulation results suggest that the Monte-Carlo version of 
\code{\link{gvtest}} statistic is more accurate
and powerful than its competitors proposed by Box and Pierce (1970), Ljung and Box (1978),
and Pena and Rodriguez (2002,2006) in the univariate time series and
Hosking (1980) and Li and McLeod (1981) in the multivariate time series.

For faster calculations, users with multi-core PC may take advantage of implementing of the 
parallel computing techniques in \code{\link{portest}} function after loading 
the package \code{snow}.
The default argument in \code{\link{portest}} function, \code{nslaves=1}, 
implements the Monte-Carlo test on PC with only one CPU. 
When the argument \code{nslaves>1} is used, 
the package \code{snow} is used provided that the 
default argument \code{MonteCarlo=TRUE} is selected.
When \code{MonteCarlo=FALSE} is selected, the test statistic 
selected from the argument \code{test} is implemented
based on the asymptotic approximation distribution.
The default test statistic is the generalized variance test, \code{\link{gvtest}}.

The \code{snow} \code{R} package is available on \code{CRAN} with more details in the 
link \url{http://www.stat.uiowa.edu/~luke/R/cluster/cluster.html}.

\subsection{Test for usual residuals and \code{ARCH} effects}{

By setting the argument \code{SquaredQ = TRUE} in \code{\link{portest}} function,
the portmanteau test based on the Monte-Carlo or the asymptotic distribution method
(depending on the choice of the argument \code{MonteCarlo})
for \code{ARCH} effects will be implemented on the squared residuals.
Otherwise, the portmanteau test will be applied on the usual residuals when the 
default argument \code{SquaredQ = FALSE} is selected.
}

\subsection{Monte-Carlo test for residuals with infinite variances}{

The argument \code{InfiniteVarianceQ = TRUE} in \code{\link{portest}} function is used only with
Monte-Carlo techniques. 
By selecting this argument, the Monte-Carlo diagnostic test on residuals with
infinite variances effects is implemented.
}

\subsection{Test for fractional Gaussian noise, \code{FGN}, effects}{

After fitting \code{FGN} model using the function \code{FitFGN} in the \code{FGN} \code{R} package,
the output object has a class \code{"FitFGN"}.
By substituting this object as a first entry in the \code{\link{portest}} function, the portmanteau
test based on the Monte-Carlo or the asymptotic distribution method (depending on the choice of \code{MonteCarlo})
for \code{FGN} model will be implemented.
}
}

\section{Simulate data from nonseasonal \code{ARIMA(p,d,q)} or \code{VARIMA(p,d,q)}}{
The function \code{\link{varima.sim}} in this package is useful for simulating data
from nonseasonal \code{ARIMA} or \code{VARIMA} of order \eqn{(p,d,q)} 
with or without deterministic terms (drift and trend), where
innovations generated from stable distribution with finite or infinite variances.
The argument \eqn{d} must entered as a nonnegative integer in the \code{ARIMA} case, whereas 
it must entered as a vector of \eqn{k} components \eqn{d_1,...,d_k} in the \code{VARIMA} case.
\eqn{d_i} represents the difference lag need to be applied on series i.
The default argument of \code{\link{varima.sim}} function, \code{StableParameters=NA}, 
is used for innovations with finite variances. 
For innovations from stable distribution with stable parameters \code{ALPHA, BETA, GAMMA,} 
and \code{DELTA}, the argument \code{StableParameters} equals to these parameters. 
}

\examples{
\dontrun{
#########################################################################
## Monte-Carlo test for randomness On PC with Single CPU ##
#########################################################################
data("DEXCAUS")
returns <- log(DEXCAUS[-1]/DEXCAUS[-length(DEXCAUS)])
portest(returns)
#########################################################################
## Monte-Carlo test for usual residuals using "snow" ##
#########################################################################
ans <- arima(Nile, order=c(1,0,1))
portest(ans, nslaves = 8)
##
IBMSP500 <- monthibmspln
Fit3 <-  ar.ols(IBMSP500, aic=TRUE,intercept=F, order.max=3)
portest(Fit3,test="gvtest",nslaves=8)
portest(Fit3,test="Hosking",nslaves=8)
## Use VAR function in the R package "vars"
##
library(vars)
Fit3 <- VAR(IBMSP500,p=3)
portest(Fit3,test="gvtest",nslaves=8)
#########################################################################
## Monte-Carlo test for squared residuals: ARCH effects using "snow" ##
#########################################################################
data("monthintel")
returns <- as.ts(monthintel)
lags <- c(5, 10, 40)
portest(returns, lags = lags, nslaves = 8, SquaredQ = TRUE)
#########################################################################
## Monte-Carlo test on residuals with infinite variances using "snow" ##
#########################################################################
data("CRSP")
CRSP.AR5<- arima(CRSP, c(5, 0, 0))
NREP<- 1000
lags=c(10, 20, 30)
portest(CRSP.AR5,lags=lags,nslaves=8,NREP=NREP,InfiniteVarianceQ=TRUE)      
#########################################################################
## Monte-Carlo test for FGN using "snow" ##
#########################################################################
library("FGN")
data("NileMin")
NILE.FGN <- FitFGN(NileMin)
lags <- c(5, 10, 20)
portest(NILE.FGN, lags = lags, nslaves = 8) 
#########################################################################
## Simulate a VARMA process with innovations that have finite variance ##
#########################################################################
## Example 1
## Data from VARMA (1,1) with mean zero vector without drift or trend.
#
n <- 300
phi <-  array(c(0.5, 0.4, 0.1, 0.5),dim = c(2, 2, 1))
theta <-  array(c(0, 0.25, 0, 0),dim = c(2, 2, 1))
d <- NA
sigma <- matrix(c(1, 0.71, 0.71, 2), 2, 2)
Z <- varima.sim(phi, theta, d, sigma, n)
plot(Z)
#
## Example 2
## Mean = c(10,20), drift equation a + b * t, a = c(2,3), and b = c(0,0.06).
#
n <- 300  
phi <-  array(c(0.5, 0.4, 0.1, 0.5),dim = c(2, 2, 1))
theta <-  array(c(0, 0.25, 0, 0),dim = c(2, 2, 1))
d <- c(0,0)
sigma <- matrix(c(1, 0.71, 0.71, 2), 2, 2)
constant <- c(2,3)
trend <- c(0,0.06)
demean <- c(10,20)
Z <- varima.sim(phi, theta, d, sigma, n, constant, trend, demean)
plot(Z)
#########################################################################
##  Simulate from VARIMA process with innovations of infinite variance ##
#########################################################################
n <- 500
phi <-  array(c(0.5, 0.4, 0.1, 0.5),dim = c(2, 2, 1))
theta <-  array(c(0, 0.25, 0, 0),dim = c(2, 2, 1))
d <- c(1,0)
sigma <- matrix(c(1, 0.71, 0.71, 2), 2, 2)
ALPHA <- c(1.3, 1.6)
BETA <- c(0, 0.2)
GAMMA <- c(1, 1)
DELTA <- c(0, 0.2)
StableParameters <- c(ALPHA, BETA, GAMMA, DELTA)
varima.sim(phi, theta, d, sigma, n, StableParameters = StableParameters)
}    
}

\author{
Author: Esam Mahdi and A. Ian McLeod.

Maintainer: A. Ian McLeod <aim@stats.uwo.ca>
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

Mahdi, E. and McLeod, A.I. (2011). "Improved Multivariate Portmanteau Diagnostic Test". Submitted.

McCulloch, J. H. (1986). "Simple Consistent Estimator of Stable Distribution Parameters".
Commun. Statist.--Simula., 15(4), 1109-1136.

McLeod A.I, Li W.K (1983). "Distribution of the Residual Autocorrelation in Multivariate
ARMA Time Series Models". Journal of Time Series Analysis, 4, 269-273.

McLeod, A.I., Yu, Hao, Krougly, Zinovi L. (2007). "Algorithms for Linear Time Series Analysis".
Journal of Statistical Software.

Pena, D. and Rodriguez, J. (2006). "The log of the determinant of the autocorrelation
matrix for testing goodness of fit in time series". Journal of Statistical Planning and
Inference, 136, 2706-2718.
}

\keyword{ PACKAGE }
\keyword{htest}
\keyword{ ts}



