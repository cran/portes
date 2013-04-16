\name{varima.sim}
\alias{varima.sim}

\title{Simulate Data From Nonseasonal ARIMA(p,d,q) or VARIMA(p,d,q) Models}

\description{
 Simulate time series from Integrated AutoRegressive Moving Average, \code{ARIMA(p,d,q)}, or 
 Vector Integrated AutoRegressive Moving Average, \code{VARIMA(p,d,q)}, 
 where \code{d} is a nonnegative difference integer in the \code{ARIMA} case 
 and it is a vector of \eqn{k} differenced components
 \eqn{d_1,...,d_k} in the \code{VARIMA} case.
 The simulated process may have a deterministic terms, drift constant and time trend, with non-zero mean.
 The innovations may have finite or infinite variance.
}

\usage{
varima.sim(phi=NULL,theta=NULL,d=NA,sigma,n,constant= NA,trend=NA, 
         demean=NA,innov=NULL,innov.dist=c("normal","t","stable"),  
         df=1,StableParameters=NA,Trunc.Series=NA)
}

\arguments{
  \item{phi}{a numeric or an array of \code{AR} or an array of \code{VAR} parameters with order \eqn{p}.}
  \item{theta}{a numeric or an array of \code{MA} or an array of \code{VMA} parameters with order \eqn{q}.}
  \item{d}{an integer or a vector representing the order of the difference.}
  \item{sigma}{variance of white noise series. It must be entered as matrix in case of bivariate or multivariate time series.}
  \item{n}{length of the series. }
  \item{constant}{a numeric vector represents the intercept in the deterministic equation.}
  \item{trend}{a numeric vector represents the slop in the deterministic equation.}
  \item{demean}{a numeric vector represents the mean of the series.}
  \item{innov}{univariate or multivariate innovation series.}
  \item{innov.dist}{innovation distribution from univariate or multivariate \code{normal}, \code{t}, or \code{stable}.
   Default is normal when \code{innov = NULL} is selected.}
  \item{df}{degrees of freedom needed for generating innovations with univariate or multivariate t-distribution.}
  \item{StableParameters}{the four parameters, \code{ALPHA, BETA, GAMMA,} and \code{DELTA}, as 
       described in the function \code{\link{rStable}}. 
       This argument is needed when both arguments \code{innov.dist = "stable"} and \code{innov = NULL} are selected. 
       These parameters are used to generate data from \code{ARIMA} or \code{VARIMA} 
       with innovations from stable distribution with infinite variance.}
  \item{Trunc.Series}{truncation lag is used to truncate the infinite \code{MA} or \code{VMA} Process. 
   IF it is \code{NA}, then \code{Trunc.Series} = min(100,\eqn{n/3}). }
}

\details{
This function is used to simulate a nonseasonal \code{ARIMA} or \code{VARIMA} model of order \eqn{(p,d,q)} 
 \deqn{\Phi(B)D(B)(Z_{t}-\mu) = a + b \times t + \Theta(B)e_{t},}
where \eqn{a, b,} and \eqn{\mu} correspond to the arguments \code{constant}, \code{trend}, and \code{demean} respectively.
The univariate or multivariate series \eqn{e_{t}} represents the innovations series given from the argument \code{innov}.
If \code{innov = NULL} then \eqn{e_{t}} will be generated from   
a univariate or multivariate normal distribution, t-distribution, or stable 
distribution of infinite variance specified from the argument \code{innov.dist}. 
\eqn{\Phi(B)} and \eqn{\Theta(B)} are the \code{VAR} \code{VMA} coefficient matrices respectively 
and \eqn{B} is the backshift time operator.
\eqn{D(B)=diag[(1-B)^{d_{1}},\ldots,(1-B)^{d_{k}}]} is a diagonal matrix. 
This states that each individual series \eqn{Z_{i}, i=1,...,k} is differenced \eqn{d_{i}} times
to reduce to a stationary \code{VARMA(p,q)} series.
}

\value{
Simulated data from \code{ARIMA(p,d,q)} or \code{VARIMA(p,d,q)} process 
that may have a drift and deterministic time trend terms.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Hipel, K.W. and McLeod, A.I. (2005). "Time Series Modelling of Water Resources and Environmental Systems".

Reinsel, G. C. (1997). "Elements of Multivariate Time Series Analysis". Springer-Verlag, 2nd edition.
}

\seealso{
\code{\link[stats]{arima.sim}}, \code{\link{vma.sim}}, \code{\link{ImpulseVMA}}, \code{\link{InvertQ}}, 
\code{\link{fitstable}}, \code{\link{rStable}}
}

\examples{
########################################################################
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
#to save time the other examples are not run
\dontrun{
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
# Simulate a bivariate VARIMA(2,d,1) process with length 300, where d=(1,2). 
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
\keyword{ ts}

