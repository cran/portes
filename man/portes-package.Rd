\name{portes-package}
\alias{portes-package}
\docType{package}

\title{
Portmanteau Tests for ARMA, VARMA, ARCH, and FGN Models
}

\description{
This package contains a set of portmanteau diagnostic checks for univariate and multivariate time series,
based on the chi-square asymptotic distribution and the parametric bootstrap significance level.
It can be used for simulating data from \code{ARMA} or \code{VARMA} Models. 
}

\details{
\tabular{ll}{
Package: \tab portes\cr
Type: \tab Package\cr
Version: \tab 1.04\cr
Date: \tab 2010-09-08\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
Depends: \tab FitAR, akima \cr
Suggests: \tab FGN, snow, rlecuyer \cr
Recommends: \tab Rmpi \cr
Classification/ACM: G.3, G.4, I.5.1\cr
Classification/MSC: 62M10, 91B84\cr
License: GPL (>= 2)\cr
}
}

\section{Parametric Bootstrap Portmanteau Test}{
The main function in this package, \code{\link{portes}}, 
works with univariate and multivariate time series.
It implements the parametric bootstrap versions of the generalized variance test as described by 
Lin and McLeod (2006), Mahdi and McLeod (2010), and the modified test given by Li and McLeod (1981)
(See the argument \code{statistic=c("GVStat","LiMcLeod")} in \code{\link{portes}} function).
The generalized variance portmanteau test is implemented in this package with the name \code{\link{GVStat}},
while the other one has the name \code{\link{LiMcLeod}}.
Both tests, \code{\link{GVStat}} and \code{\link{LiMcLeod}} can be used for univariate and multivariate time series.   
For faster calculations, users that have PC's with multiple CPU's may take advantage of implementing the 
parallel computing techniques in \code{\link{portes}} function after downloading either \code{Rmpi} or \code{snow}.
The default argument in \code{\link{portes}} function, \code{SpawnSlaves=c("default","Rmpi","snow")},
implements the parametric bootstrap test on PC's with only one CPU's.
The \code{snow} \code{R} package is available on \code{CRAN} with more details in the 
link \url{http://www.stat.uiowa.edu/~luke/R/cluster/cluster.html}.
Instructions to install and run \code{Rmpi} package under Windows is given in the link
\url{http://www.stats.uwo.ca/faculty/yu/Rmpi/windows.htm}.
More information about parallel computing with \code{R} is provided by 
\url{http://cran.r-project.org/web/views/HighPerformanceComputing.html}.

\subsection{parametric bootstrap test for usual residuals}{

The default arguments in the function \code{\link{portes}} corresponds to  
parametric bootstrap test for the usual residuals. 
}

\subsection{parametric bootstrap test for squared residuals, \code{ARCH} effects}{

By setting the argument \code{SquaredQ = TRUE} in \code{\link{portes}} function,
the parametric bootstrap test for \code{ARCH} effects will be implemented.
}

\subsection{parametric bootstrap test for residuals with infinite variances}{

By setting the argument \code{InfiniteVarianceQ = TRUE} in \code{\link{portes}} function,
the parametric bootstrap test for residuals with infinite variance will be implemented.
}

\subsection{parametric bootstrap test for \code{FGN}}{

After fitting \code{FGN} model using the function \code{FitFGN} in the \code{FGN} \code{R} package,
the output will be an object with a class \code{"FitFGN"}.
By substituting this object as a first entry in the \code{\link{portes}} function, the
parametric bootstrap test for \code{FGN} model will be implemented.
}
}

\section{Simulate Data From \code{ARMA} or \code{VARMA}}{
The function \code{\link{simvarma}} in this package is very useful for simulating data
from \code{ARMA} or \code{VARMA} with innovations follow Gaussian or stable Paretian distribution.
The default argument of \code{\link{simvarma}} function, \code{StableParameters=NA}, is used for Gaussian innovations. 
For innovations from stable Paretian distribution with stable parameters \code{ALPHA, BETA, GAMMA,} 
and \code{DELTA}, the argument \code{StableParameters} is setting to equal to these parameters values. 
}

\examples{
\dontrun{
#########################################################################
## parametric bootstrap test for usual residuals ##
#########################################################################
data("DEXCAUS")
returns <- log(DEXCAUS[-1]/DEXCAUS[-length(DEXCAUS)])
portes(returns)

data("house")
HouseData <- scale(house, scale = FALSE)
HOUSE.VAR11<- ar.ols(HouseData,aic=FALSE,intercept=FALSE,order.max=11)
res <- ts((HOUSE.VAR11$resid)[-(1:11), ])
portes(res, SpawnSlaves="snow")
#########################################################################
## parametric bootstrap test for squared residuals, ARCH effects ##
#########################################################################
data("monthintel")
returns <- as.ts(monthintel)
lags <- c(5, 10, 40)
portes(returns,lags,"GVStat",SpawnSlaves="Rmpi",nslaves=8,SquaredQ=TRUE)
#########################################################################
## parametric bootstrap test for residuals with infinite variances ##
#########################################################################
data("CRSP")
CRSP.AR5<- arima(CRSP, c(5, 0, 0))
NREP<- 1000
lags=c(10, 20, 30)
portes(CRSP.AR5, lags = lags, statistic = "GVStat", SpawnSlaves = "Rmpi",
        nslaves = 8, NREP = NREP, InfiniteVarianceQ = TRUE)      
#########################################################################
## parametric bootstrap test for FGN ##
#########################################################################
library("FGN")
data("NileMin")
NILE.FGN <- FitFGN(NileMin)
lags <- c(5, 10, 20)
portes(NILE.FGN, lags, "GVStat", SpawnSlaves = "Rmpi", nslaves = 8) 
#########################################################################
## Simulate a VARMA process with innovations that have finite variance ##
#########################################################################
n <- 300
phi <-  array(c(0.5, 0.4, 0.1, 0.5),dim = c(2, 2, 1))
theta <-  array(c(0, 0.25, 0, 0),dim = c(2, 2, 1))
intercept <- c(1, 2)
sigma <- matrix(c(1, 0.71, 0.71, 2), 2, 2)
simvarma(phi, theta, sigma, intercept, n)
#########################################################################
## Simulate a VARMA process with innovations that have infinite variance#
#########################################################################
n <- 300
phi <-  array(c(0.5, 0.4, 0.1, 0.5),dim = c(2, 2, 1))
theta <-  array(c(0, 0.25, 0, 0),dim = c(2, 2, 1))
intercept <- c(1, 2)
sigma <- matrix(c(1, 0.71, 0.71, 2), 2, 2)
ALPHA <- c(1.3, 1.6)
BETA <- c(0, 0.2)
GAMMA <- c(1, 1)
DELTA <- c(0, 0.2)
StableParameters <- c(ALPHA, BETA, GAMMA, DELTA)
simvarma(phi, theta, sigma, intercept, n, StableParameters)
}    
}

\author{
Author: Esam Mahdi & A. Ian McLeod.

Maintainer: A. Ian McLeod <aim@stats.uwo.ca>
}


\references{
Li, W. K. and McLeod, A. I. (1981). "Distribution of The Residual Autocorrelations in 
Multivariate ARMA Time Series Models". Journal of The Royal Statistical Society, Series B, 43, 231-239.

Lin, J.-W. and McLeod, A.I. (2006). "Improved Generalized Variance Portmanteau Test".
Computational Statistics and Data Analysis 51, 1731-1738.

Lin, J.-W. and McLeod, A.I. (2008). "Portmanteau Tests for ARMA Models with Infinite Variance".  
Journal of Time Series Analysis, 29, 600-617.

Mahdi, E. and McLeod, A.I. (2010). "Improved Multivariate Portmanteau Diagnostic Test". Submitted.

McCulloch, J. H. (1986). "Simple Consistent Estimator of Stable Distribution Parameters".
Commun. Statist.--Simula., 15(4), 1109-1136.

McLeod AI, Li WK (1983). "Distribution of the Residual Autocorrelation in Multivariate
ARMA Time Series Models". Journal of Time Series Analysis, 4, 269-273.

McLeod, A.I., Yu, Hao, Krougly, Zinovi L. (2007). "Algorithms for Linear Time Series Analysis".
Journal of Statistical Software.
}

\keyword{ PACKAGE }
\keyword{htest}
\keyword{ ts}



