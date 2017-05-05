\name{vma.sim}
\alias{vma.sim}

\title{Compute The Vector of Moving Average Model (VMA)}

\description{
This utility function is useful to use in the function \code{\link{varima.sim}} 
and may used to compute the coefficients of moving-average or vector moving-average.
}

\usage{
vma.sim(psi, a)
}

\arguments{
  \item{psi}{the impulse coefficients.}
  \item{a}{innovations}
}

\value{
Vector of length \eqn{n} (in the univariate case), or \eqn{n} matrices (in the multivariate case), 
where \eqn{n} = length(\eqn{a})-length(\eqn{\Psi}) and \eqn{n\times k} is the dimension of the series.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Hannan, E.J. (1970). "Multiple Time Series". New York: Wiley.

Hipel, K.W. and McLeod, A.I. (2005). "Time Series Modelling of Water Resources and Environmental Systems".
}

\seealso{
\code{\link[stats]{convolve}}, \code{\link{varima.sim}}, \code{\link[stats]{arima.sim}}, \code{\link{ImpulseVMA}}, 
\code{\link{InvertQ}}, \code{\link{fitstable}}
 }

\examples{
k <- 2
n <- 300
trunc.lag <-  50   
phi <-  array(c(0.5,0.4,0.1,0.5),dim=c(k,k,1))
theta <-  array(c(0,0.25,0,0),dim=c(k,k,1))
sigma <- matrix(c(1,0.71,0.71,2),k,k)
p <- ifelse(is.null(phi),0,dim(phi)[3])
q <- ifelse(is.null(theta),0,dim(theta)[3])
r <- max(p, q)
d <- trunc.lag + r
psi <- ImpulseVMA(phi = phi, theta = theta, trunc.lag = trunc.lag)
a <- t(crossprod(chol(sigma),matrix(rnorm(k*d),ncol=d)))
vma.sim(psi = psi, a = a)
}

\keyword{ ts}

