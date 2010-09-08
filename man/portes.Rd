\name{portes}
\alias{portes}

\title{Parametric Bootstrap Portmanteau Tests}

\description{
 Univariate or multivariate \code{\link{GVStat}} and \code{\link{LiMcLeod}} portmanteau 
tests using parametric bootstrap techniques.
}

\usage{
  portes(obj, lags=seq(5,30,5), statistic=c("GVStat","LiMcLeod"), 
  order=0, SpawnSlaves=c("default","Rmpi","snow"), nslaves=2, 
  NREP=1000, InfiniteVarianceQ=FALSE,SquaredQ=FALSE, Set.Seed=TRUE)
}

\arguments{
  \item{obj}{if \code{obj} is an object of class \code{"ar"}, \code{"Arima"}, \code{"arima0"},
   \code{"FitAR"}, or \code{"FitFGN"} then a portmanteau goodness-of-fit test is done on the 
    fitted model. 
    Otherwise, for \code{obj} with class \code{"ts"}, \code{"numeric"}, \code{"matrix"}, 
    or \code{("mts" "ts")}, a test of randomness is done.}
  \item{lags}{vector of lag values is used for portmanteau test.}
  \item{statistic}{character string giving the name of the portmanteau test, default is "\code{\link{GVStat}}".}
  \item{order}{\eqn{p+q}, where (\eqn{p},\eqn{q}) is the \code{ARMA} or the \code{VARMA} order.}
  \item{SpawnSlaves}{\code{"default"} when a computer with single CPU is used, 
       \code{"Rmpi"} or \code{"snow"} when a computer with MultiCore CPUs is used.}
  \item{nslaves}{number of slaves needed to use in parallel calculations.}
  \item{NREP}{number of replications needed for parametric bootstrap test.}
  \item{InfiniteVarianceQ}{\code{FALSE}, assumes innovations with finite variance from Gaussian distribution. 
    Otherwise, innovations follow stable Paretian distribution with infinite variance.}
  \item{SquaredQ}{\code{TRUE}, apply the test to the squared values. 
    This is needed to check the Autoregressive Conditional Heteroscedastic, \code{ARCH}, effects. 
    When \code{SquaredQ = FALSE}, apply for the usual residuals. }
  \item{Set.Seed}{initialize random seed.}
}

\details{
The generalized variance portmanteau statistic, \code{\link{GVStat}}, and the modified portmanteau statistic, 
\code{\link{LiMcLeod}}, are implemented based on parametric bootstrap techniques. 
The null hypothesis assuming that the fitted model is an adequate model and the residuals behave like white noise series. 
The \code{Rmpi} package must be installed if one decide to choose the argument \code{SpawnSlaves="Rmpi"} and
the \code{snow} package must be installed if one decide to choose the argument \code{SpawnSlaves="snow"}.
This function can be used for testing \code{ARCH} effects, and for checking  
adequacy in the fitted \code{ARMA}, \code{VAR}, and Fractional Gaussian Noise, \code{FGN}, models.
Also, it can be used to check for randomness of residuals.
}

\value{
vector of p-values for different lag values.
}

\author{ Esam Mahdi & A.I. McLeod.}

\references{
Li, W. K. and McLeod, A. I. (1981). "Distribution of The Residual Autocorrelations in 
Multivariate ARMA Time Series Models". Journal of The Royal Statistical Society, Series B, 43, 231-239.

Lin, J.-W. and McLeod, A.I. (2006). "Improved Generalized Variance Portmanteau Test".
Computational Statistics and Data Analysis 51, 1731-1738.

Lin, J.-W. and McLeod, A.I. (2008). "Portmanteau Tests for ARMA Models with Infinite Variance".  
Journal of Time Series Analysis, 29, 600-617.

Mahdi, E. and McLeod, A.I. (2010). "Improved Multivariate Portmanteau Diagnostic Test". Submitted.

McLeod, A.I., Yu, Hao, Krougly, Zinovi L. (2007). "Algorithms for Linear Time Series Analysis".
Journal of Statistical Software.
}

\seealso{
 \code{\link{simvarma}}, \code{\link{ar}}, \code{\link{arima}}, \code{\link{arima0}}, \code{FitAR}, \code{FitFGN},
\code{\link{LiMcLeod}}, \code{\link{GVStat}}, \code{\link{FitStable}}
 }

\examples{
#Example 1
#Simple example using portes.
# 100 replications takes about 6 seconds on single cpu 2.67 GHz
ans<-arima(Nile, order=c(1,0,1))
portes(ans, NREP=100)
#
#Example 2
#Simulate a bivariate VARMA(1,0) process with length 300. 
k <- 2
n <- 300
Trunc.Series <-  50   
phi <-  array(c(0.5,0.4,0.1,0.5),dim=c(k,k,1))
theta <-  array(c(0,0.0,0,0),dim=c(k,k,1))
intercept <- c(10,12)
sigma <- matrix(c(1,0.71,0.71,2),k,k)
z<-simvarma(phi,theta=theta,sigma,intercept,n)
ans<-ar(z)
portes(ans, NREP=25) #NREP=25 is not enough for accurate p-values!
#
#Examples using snow or Rmpi. Takes a few minutes.
\dontrun{
##############################################################
## The following script checks the residuals for randomness 
## Model VAR(11) fitted to house data after subtracting the mean. 
##
##
HouseData <- scale(house, scale = FALSE)
HOUSE.VAR11<-ar.ols(HouseData,aic=FALSE,intercept=FALSE,order.max=11)
res <- ts((HOUSE.VAR11$resid)[-(1:11), ])
#
## GVStat statistic on PC with one CPU (61.68 seconds)
portes(res) 
# 
## GVStat statistic on PC with 8 CPU's using "Rmpi" (10.44 seconds)
portes(res,statistic="GVStat",SpawnSlaves="Rmpi", nslaves=8)
#
## GVStat statistic on PC with 8 CPU's using "snow" (10.81 seconds)
portes(res,statistic="GVStat",SpawnSlaves="snow", nslaves=8)
##############################################################
## The following script generates Table 2 in Lin and McLeod (2008).
## It takes 1.96 minutes on PC with 8 CPU's using "Rmpi"
##
CRSP.AR5<- arima(CRSP, c(5,0,0))
NREP<- 1000
lags=c(10,20,30)
tb<-matrix(numeric(6*length(lags)),nrow=6)
tb[1,]<-portes(CRSP.AR5,lags=lags,statistic="GVStat",SpawnSlaves="Rmpi",
        nslaves=8, NREP=NREP,InfiniteVarianceQ=TRUE)[,2]
tb[2,]<-portes(CRSP.AR5,lags=lags,statistic="LiMcLeod",SpawnSlaves="Rmpi", 
        nslaves=8, NREP=NREP,InfiniteVarianceQ=TRUE)[,2]
tb[3,]<-portes(CRSP.AR5,lags=lags,statistic="GVStat",SpawnSlaves="Rmpi",
        nslaves=8, NREP=NREP)[,2]
tb[4,]<-portes(CRSP.AR5,lags=lags,statistic="LiMcLeod",SpawnSlaves="Rmpi",
        nslaves=8, NREP=NREP)[,2]
tb[5,]<-GVStat(CRSP.AR5,lags,5)[,4]
tb[6,]<-LiMcLeod(CRSP.AR5,lags,5)[,4]
rn <- c("GVStat-MC-Stable","LiMcLeod-MC-Stable", "GVStat-MC-Normal",
        "LiMcLeod-MC-Normal", "GVStat-asymptotic", "LiMcLeod-asymptotic")
dimnames(tb) <- list(rn, lags)
tb
#
#Output
#                            10          20          30
#GVStat-MC-Stable    0.06393606 0.034965035 0.021978022
#LiMcLeod-MC-Stable  0.05494505 0.019980020 0.022977023
#GVStat-MC-Normal    0.06293706 0.014985015 0.001998002
#LiMcLeod-MC-Normal  0.04995005 0.001998002 0.001998002
#GVStat-asymptotic   0.29111084 0.057389935 0.007760695
#LiMcLeod-asymptotic 0.04858500 0.002316690 0.002643823
##############################################################
## Fit AR(2) to Nile river data on PC with 8 CPU's
##
##
NILE.AR2 <- FitAR(Nile, 2)
lags <- c(5, 10, 20)
##				
## GVStat test on fitted model using "Rmpi" (5.9 seconds)
portes(res, lags, "GVStat",SpawnSlaves="Rmpi", nslaves=8) 
#
## LiMcLeod test on residuals using "snow" (5.95 seconds)
res <- NILE.AR2$res
portes(NILE.AR2, lags, "LiMcLeod",SpawnSlaves="snow", nslaves=8) 
##############################################################
## Test monthly log stock returns of Intel data for ARCH effects
## GVStat statistic on PC with 8 CPU's using "Rmpi" 
## It takes 5.93 seconds to get the results
##
data(monthintel)
returns <- as.ts(monthintel)
lags <- c(5, 10, 40)
portes(returns,lags,"GVStat",SpawnSlaves="Rmpi",nslaves=8,SquaredQ=TRUE) 
##############################################################
## Fit Fractional Gaussian Noise, FGN, to NileMin data in FGN package
## GVStat statistic on PC with 8 CPU's using "Rmpi"
##
library(FGN)
data(NileMin)
NILE.FGN <- FitFGN(NileMin)
lags <- c(5, 10, 20)
##				
## GVStat statistic on fitted model (55 seconds)
portes(NILE.FGN, lags, "GVStat",SpawnSlaves="Rmpi", nslaves=8) 
#
## GVStat statistic on residuals (6 seconds)
res <- NILE.FGN$res
portes(res, lags, "GVStat",SpawnSlaves="Rmpi", nslaves=8)
}
}

\keyword{htest}
\keyword{ ts}


