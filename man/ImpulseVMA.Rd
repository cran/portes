\name{ImpulseVMA}
\alias{ImpulseVMA}

\title{The Impulse Response Function in the Infinite MA or VMA Representation}

\description{
The impulse coefficients are computed.
}

\usage{
ImpulseVMA(phi=NULL,theta=NULL,trunc.lag=NULL)
}

\arguments{
  \item{phi}{a numeric or an array of \code{AR} or an array of \code{VAR} parameters with order \eqn{p}.}
  \item{theta}{a numeric or an array of \code{MA} or an array of \code{VMA} parameters with order \eqn{q}.}
  \item{trunc.lag}{truncation lag is used to truncate the infinite \code{MA} or \code{VMA} Process. 
  IF it is \code{NULL}, then the default \code{trunc.lag} = \eqn{p+q}.}
}

\value{
The impulse response coefficients of order \code{trunc.lag+1} obtained by 
converting the \code{ARMA}\eqn{(p,q)} or \code{VARMA}\eqn{(p,q)} 
process to infinite \code{MA} or \code{VMA} process, respectively. 
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Lutkepohl, H. (2005). "New introduction to multiple time series analysis". 
Springer-Verlag, New York.

Reinsel, G. C. (1997). "Elements of Multivariate Time Series Analysis". 
Springer-Verlag, 2nd edition.
}

\seealso{ \code{\link[stats]{ARMAtoMA}}, \code{\link{varima.sim}}, \code{\link{vma.sim}}, 
\code{\link{InvertQ}}, \code{\link[FitAR]{InvertibleQ}}
}

\examples{
#####################################################################
### Impulse response coefficients from AR(1,1) to infinite MA process. 
### The infinite process is truncated at lag 20
###
k <- 1
trunc.lag <- 20
phi <- 0.7
theta <- array(-0.9,dim=c(k,k,1))
ImpulseVMA(phi,theta,trunc.lag)
#####################################################################
### Impulse response coefficients from VAR(2) to infinite VMA process
### The infinite process is truncated at default lag value = p+q
###
k <- 2
phi <- array(c(0.5,0.4,0.1,0.5,0,0.3,0,0),dim=c(k,k,2))
theta <- NULL
ImpulseVMA(phi,theta)
#####################################################################
### Impulse response coefficients from VARMA(2,1) to infinite VMA process
### The infinite process is truncated at lag 50
###
k <- 2
phi <- array(c(0.5,0.4,0.1,0.5,0,0.25,0,0),dim=c(k,k,2))
theta <- array(c(0.6,0,0.2,0.3),dim=c(k,k,1))
ImpulseVMA(phi,theta,trunc.lag=50)
}

\keyword{ ts}

