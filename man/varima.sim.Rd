\name{varima.sim}
\alias{varima.sim}

\title{Simulate Data From Seasonal/Nonseasonal ARIMA(p,d,q)*(ps,ds,qs)_s or VARIMA(p,d,q)*(ps,ds,qs)_s Models}

\description{
 Simulate time series from AutoRegressive Integrated Moving Average, \code{ARIMA(p,d,q)}, or 
Vector Integrated AutoRegressive Moving Average, \code{VARIMA(p,d,q)}, 
 where \code{d} is a nonnegative difference integer in the \code{ARIMA} case 
 and it is a vector of \eqn{k} differenced components
 \eqn{d_1,...,d_k} in the \code{VARIMA} case.
 In general, this function can be implemented in simulating univariate or multivariate 
 Seasonal AutoRegressive Integrated Moving Average, \code{SARIMA(p,d,q)*(ps,ds,qs)_s} 
 and \code{SVARIMA(p,d,q)*(ps,ds,qs)_s}, where \code{ps} and \code{qs} are 
 the orders of the seasonal univariate/multivariate AutoRegressive and Moving Average components respectively.
 \code{ds} is a nonnegative difference integer in the \code{SARIMA} case 
 and it is a vector of \eqn{k} differenced components \eqn{ds_1,...,ds_k} in the \code{SVARIMA} case, 
 whereas \code{s} is the seasonal period. 
 The simulated process may have a deterministic terms, drift constant and time trend, with non-zero mean.
 The innovations may have finite or infinite variance.
}

\usage{
varima.sim(model=list(ar=NULL,ma=NULL,d=NULL,sar=NULL,sma=NULL,D=NULL,period=NULL),
          n,k=1,constant=NA,trend=NA,demean=NA,innov=NULL, 
          innov.dist=c("Gaussian","t","bootstrap"),...)
}

\arguments{
  \item{model}{a list with univariate/multivariate component \code{ar} and/or \code{ma} and/or \code{sar}
   and/or \code{sma} giving the univariate/multivariate \code{AR} and/or \code{MA} and/or \code{SAR} 
   and/or \code{SMA} coefficients respectively. 
  \code{period} specifies the seasonal period. 
  For seasonality, default is \code{NULL} indicates that \code{period =12}. 
  \code{d} and \code{D} are integer or vector representing the order of the usual and seasonal difference.
  An empty list gives an \code{ARIMA(0, 0, 0)*(0,0,0)_null} model, that is white noise.}
  \item{n}{length of the series.}
  \item{k}{ number of simulated series. For example, \code{k=1} is used for univariate series and \code{k=2} is used for bivariate series.} 
  \item{constant}{a numeric vector represents the intercept in the deterministic equation.}
  \item{trend}{a numeric vector represents the slop in the deterministic equation.}
  \item{demean}{a numeric vector represents the mean of the series.}
  \item{innov}{a vector of univariate or multivariate innovation series. 
  This may used as an initial series to genrate innovations with \code{innov.dist = "bootstrap"}.
  This argument is irrelevant with the other selections of \code{innov.dist}.}
  \item{innov.dist}{distribution to generate univariate or multivariate innovation process. 
   This could be \code{Gaussian}, \code{t}, or \code{bootstrap} using resampled errors rather than distributed errors.
   Default is \code{Gaussian}.}
  \item{...}{arguments to be passed to methods, such as \code{dft} degrees of freedom 
  needed to generate innovations with univariate/multivariate series with t-distribution innovations.
   The argument \code{trunc.lag} represents the truncation lag that is used 
   to truncate the infinite \code{MA} or \code{VMA} Process. 
   IF it is not given, then \code{trunc.lag} = min(100, \eqn{n/3}). 
   Optionally \code{sigma} is the variance of a Gaussian or t white noise series.}
}

\details{
This function is used to simulate a univariate/multivariate seasonal/nonseasonal 
\code{SARIMA} or \code{SVARIMA} model of order \eqn{(p,d,q)\times(ps,ds,qs)_s} 
 \deqn{\phi(B)\Phi(B^s)d(B)D(B^s)(Z_{t})-\mu = a + b \times t + \theta(B)\Theta(B^s)e_{t},}
where \eqn{a, b,} and \eqn{\mu} correspond to the arguments \code{constant}, \code{trend}, and \code{demean} respectively.
The univariate or multivariate series \eqn{e_{t}} represents the innovations series given from the argument \code{innov}.
If \code{innov = NULL} then \eqn{e_{t}} will be generated from   
a univariate or multivariate normal distribution or t-distribution. 
\eqn{\phi(B)} and \eqn{\theta(B)} are the \code{VAR} and the \code{VMA} coefficient matrices respectively 
and \eqn{B} is the backshift time operator.
\eqn{\Phi(B^s)} and \eqn{\Theta(B^s)} are the \code{Vector SAR} \code{Vector SMA} 
coefficient matrices respectively.
\eqn{d(B)=diag[(1-B)^{d_{1}},\ldots,(1-B)^{d_{k}}]} and
\eqn{D(B^s)=diag[(1-B^s)^{ds_{1}},\ldots,(1-B^s)^{ds_{k}}]} are diagonal matrices. 
This states that each individual series \eqn{Z_{i}, i=1,...,k} is differenced \eqn{d_{i}ds_{i}} times
to reduce to a stationary \code{Vector ARMA(p,0,q)*(ps,0,qs)_s} series.
}

\value{
Simulated data from \code{SARIMA(p,d,q)} or \code{SVARIMA(p,d,q)*(ps,ds,qs)_s} process 
that may have a drift and deterministic time trend terms.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Hipel, K.W. and McLeod, A.I. (2005). "Time Series Modelling of Water Resources and Environmental Systems".

Reinsel, G. C. (1997). "Elements of Multivariate Time Series Analysis". Springer-Verlag, 2nd edition.
}

\seealso{
\code{\link[stats]{arima.sim}}, \code{\link{vma.sim}}, \code{\link{ImpulseVMA}}, \code{\link{InvertQ}}
}

\examples{
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
ma <-  0.6
Z<-varima.sim(list(ma=ma),n=n,innov=rnorm(n),innov.dist="bootstrap")
plot(Z)
#################################################################################
# Simulate seasonal ARIMA(2,1,0)*(0,2,1)_12 process with ar=c(1.3,-0.35),      #
# ma.season = 0.8. Gaussian innovations. The series is truncated at lag 50   #
#################################################################################
set.seed(12834)
n <- 100
ar <- c(1.3, -0.35)
ma.season <- 0.8
Z<-varima.sim(list(ar=ar,d=1,sma=ma.season,D=2),n=n,trunc.lag=50)
plot(Z)
acf(Z)
#################################################################################
# Simulate seasonal ARMA(0,0,0)*(2,0,0)_4 process with intercept = 0.8          #
# ar.season = c(0.9,-0.9), period = 4, t5-distribution innovations with df = 3 #
#################################################################################
set.seed(1234)
n <- 200
ar.season <- c(0.9,-0.9)
Z<-varima.sim(list(sar=ar.season,period=4),n=n,constant=0.8,innov.dist="t",dft=3)
plot(Z)
acf(Z)
arima(Z,order=c(0,0,0),seasonal = list(order = c(2,0,0),period=4))
#################################################################################
# Simulate a bivariate white noise series from a multivariate t4-distribution   #
# Then use the nonparametric bootstrap method to generate a seasonal SVARIMA    #
# of order (0,d,0)*(0,0,1)_12 with d = c(1, 0), n= 250, k = 2, and              #
# ma.season=array(c(0.5,0.4,0.1,0.3),dim=c(k,k,1))                           #
#################################################################################
set.seed(1234)
Z1 <- varima.sim(n=250,k=2,innov.dist="t",dft=4)
ma.season=array(c(0.5,0.4,0.1,0.3),dim=c(2,2,1)) 
Z2 <- varima.sim(list(sma=ma.season,d=c(1,0)),n=250,k=2,
                 innov=Z1,innov.dist="bootstrap")
plot(Z2)
#################################################################################
# Simulate a bivariate VARIMA(2,d,1) process with length 300, where d=(1,2).    #
# ar = array(c(0.5,0.4,0.1,0.5,0,0.3,0,0),dim=c(k,k,2)),                       #
# ma = array(c(0,0.25,0,0), dim=c(k,k,1)).                                   #
# innovations are generated from multivariate normal                            #
# The process have mean zero and no deterministic terms.                        #
# The variance covariance is sigma = matrix(c(1,0.71,0.71,2),2,2).              #
# The series is truncated at default value: trunc.lag=ceiling(100/3)=34         #
#################################################################################
set.seed(1234)
k <- 2
n <- 300 
ar <-  array(c(0.5,0.4,0.1,0.5,0,0.3,0,0),dim=c(k,k,2))
ma <-  array(c(0,0.25,0,0),dim=c(k,k,1))
d <- c(1,2)
sigma <- matrix(c(1,0.71,0.71,2),k,k)
Z <- varima.sim(list(ma=ar,ma=ma,d=d),n=n,k=2,sigma=sigma)
plot(Z)
#################################################################################
# Simulate a trivariate Vector SVARMA(1,0,0)*(1,0,0)_12 process with length 300 #
# ar = array(c(0.5,0.4,0.1,0.5,0,0.3,0,0,0.1), dim=c(k,k,1)), where k =3       #
# ar.season = array(c(0,0.25,0,0.5,0.1,0.4,0,0.25,0.6), dim=c(k,k,1)).         #
# innovations are generated from multivariate normal distribution               #
# The process have mean c(10, 0, 12),                                           #
# Drift equation a + b * t, where a = c(2,1,5), and b = c(0.01,0.06,0)          #
# The series is truncated at default value: trunc.lag=ceiling(100/3)=34         #
#################################################################################
set.seed(1234)
k <- 3
n <- 300
ar <-  array(c(0.5,0.4,0.1,0.5,0,0.3,0,0,0.1),dim=c(k,k,1))
ar.season <-  array(c(0,0.25,0,0.5,0.1,0.4,0,0.25,0.6),dim=c(k,k,1))
constant <- c(2,1,5)
trend <- c(0.01,0.06,0)
demean <- c(10,0,12)
Z <- varima.sim(list(ar=ar,sar=ar.season),n=n,k=3,constant=constant,
trend=trend,demean=demean)
plot(Z)
acf(Z)
}

\keyword{ ts}

