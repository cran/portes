\name{simvarma}
\alias{simvarma}

\title{Simulate Data From ARMA(p,q) or VARMA(p,q) Models}

\description{
 Simulate time series from AutoRegressive Moving Average, \code{ARMA}\eqn{(p,q)}, or 
 Vector AutoRegressive Moving Average, \code{VARMA}\eqn{(p,q)}.
}

\usage{
simvarma(phi=NULL,theta=NULL,sigma,intercept=NA,n,
         StableParameters=NA,Trunc.Series=NA)
}

\arguments{
  \item{phi}{a numeric or an array of \code{AR} or an array of \code{VAR} parameters with order \eqn{p}.}
  \item{theta}{a numeric or an array of \code{MA} or an array of \code{VMA} parameters with order \eqn{q}.}
   \item{sigma}{variance of white noise series. It must be entered as matrix in case of bivariate or multivariate time series.}
    \item{intercept}{the mean vector of the series.}
     \item{n}{length of the series. }
  \item{StableParameters}{the four parameters, \code{ALPHA, BETA, GAMMA,} and \code{DELTA}, as 
   described in the function \code{\link{rstable}}. 
   This argument is needed to generate model with innovations assumed to have stable Paretian distribution with
    infinite variance. The default is for innovations with Gaussian distribution that have finite variance.}
     \item{Trunc.Series}{truncation lag is used to truncate the infinite \code{MA} or \code{VMA} Process. 
     IF it is \code{NA}, then \code{Trunc.Series} = min(100,\eqn{n/3}). }
}

\value{
Simulated data from \code{ARMA}\eqn{(p,q)} or \code{VARMA}\eqn{(p,q)} process with errors that have stable Paretian or Gaussian
distribution.
}

\author{ Esam Mahdi & A.I. McLeod.}

\references{
Hipel, K.W. and McLeod, A.I. (2005). "Time Series Modelling of Water Resources and Environmental Systems".

Reinsel, G. C. (1997). "Elements of Multivariate Time Series Analysis". Springer-Verlag, 2nd edition.
}

\seealso{
\code{\link[stats]{arima.sim}}, \code{\link{simvma}}, \code{\link{ImpulseVMA}}, \code{\link{InvertQ}}, 
\code{\link{FitStable}}, \code{\link{rstable}}
}

\examples{
########################################################################
# Simulate AR(2) process with phi = c(1.3, -0.35), Gaussian innovations
# The series is truncated at lag 40
Trunc.Series <- 40
n <- 200
sigma <- 1
intercept <- NA
phi <- c(1.3, -0.35)
theta <- NULL
simvarma(phi, theta, sigma, intercept, n, StableParameters=NA,Trunc.Series)
########################################################################
# Simulate a bivariate VARMA(1,1) process with length 300. 
# phi = array(c(0.5,0.4,0.1,0.5), dim=c(k,k,1)),
# theta = array(c(0,0.25,0,0), dim=c(k,k,1)).
# The innovations have normal distribution with mean c(10,12)
# The variance covariance is sigma = matrix(c(1,0.71,0.71,2),2,2).
# The series is truncated at default value: Trunc.Series=ceiling(100/3)=34 
k <- 2
n <- 300
Trunc.Series <-  50   
phi <-  array(c(0.5,0.4,0.1,0.5),dim=c(k,k,1))
theta <-  array(c(0,0.25,0,0),dim=c(k,k,1))
intercept <- c(10,12)
sigma <- matrix(c(1,0.71,0.71,2),k,k)
simvarma(phi,theta,sigma,intercept,n)
########################################################################
# Simulate univariate ARMA(2,1) process with length 500, 
# phi = c(1.3, -0.35), theta = 0.1.
# Stable innovations with: ALPHA = 1.75, BETA = 0, GAMMA = 1, DELTA = 0
n <- 500
phi <-  c(1.3, -0.35)
theta <-  0.1
intercept <- 0
sigma <-  0.7
ALPHA <- 1.75
BETA <- 0
GAMMA <- 1
DELTA <- 0
StableParameters <- c(ALPHA,BETA,GAMMA,DELTA)
simvarma(phi,theta,sigma,intercept,n,StableParameters)
########################################################################
# Simulate a bivariate VAR(1) process with length 100. 
# Stable Paretian: ALPHA=(1.3,1.6), BETA=(0,0.2), GAMMA=(1,1), DELTA=(0,0.2)
# The series is truncated at default value: Trunc.Series=min(100,200)=100 
k <- 2
n <- 600
phi <- array(c(-0.2,-0.6,0.3,1.1),dim=c(k,k,1))
theta <- NULL
sigma <- matrix(c(1,0.71,0.71,2),k,k)
intercept <- rep(0,k)
ALPHA <- c(1.3,1.6)
BETA <- c(0,0.2)
GAMMA <-c(1,1)
DELTA <-c(0,0.2)
StableParameters <- c(ALPHA,BETA,GAMMA,DELTA)
simvarma(phi,theta,sigma,intercept,n,StableParameters)
}

\keyword{ ts}

