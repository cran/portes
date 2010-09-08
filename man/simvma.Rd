\name{simvma}
\alias{simvma}

\title{Compute The Vector of Moving Average Model (VMA)}

\description{
This utility function may used to compute the coefficients of moving-average or vector moving-average.
It's also a function needed in the function \code{\link{simvarma}}. 
}

\usage{
simvma(psi, a)
}

\arguments{
  \item{psi}{the impulse coefficients.}
  \item{a}{innovations}
}

\value{
Vector of length \eqn{n} (in the univariate case), or \eqn{n} matrices (in the multivariate case), 
where \eqn{n} = length(\eqn{a})-length(\eqn{\Psi}) and \eqn{n\times k} is the dimension of the series.
}

\author{ Esam Mahdi & A.I. McLeod.}

\references{
Hannan, E.J. (1970). "Multiple Time Series". New York: Wiley.

Hipel, K.W. and McLeod, A.I. (2005). "Time Series Modelling of Water Resources and Environmental Systems".
}

\seealso{
\code{\link[stats]{convolve}}, \code{\link{simvarma}}, \code{\link[stats]{arima.sim}}, \code{\link{ImpulseVMA}}, 
\code{\link{InvertQ}}, \code{\link{FitStable}}
 }

\examples{
k <- 2
n <- 300
Trunc.Series <-  50   
phi <-  array(c(0.5,0.4,0.1,0.5),dim=c(k,k,1))
theta <-  array(c(0,0.25,0,0),dim=c(k,k,1))
intercept <- c(10,12)
sigma <- matrix(c(1,0.71,0.71,2),k,k)
p <- ifelse(is.null(phi),0,dim(phi)[3])
q <- ifelse(is.null(theta),0,dim(theta)[3])
epsilon <- t(crossprod(chol(sigma),matrix(rnorm(k*(n + q)),ncol=n + q)))
Sim.Series <- matrix(numeric(0), nrow = n, ncol = k)
r <- max(p, q)
psi <- ImpulseVMA(phi = phi, theta = theta, Trunc.Series = Trunc.Series)
epsilon <- t(crossprod(chol(sigma),matrix(rnorm(k*(Trunc.Series + r)),ncol=Trunc.Series + r)))
Sim.Series[1:r, ] <- simvma(psi = psi, a = epsilon)
Sim.Series[1:r,]
}

\keyword{ ts}

