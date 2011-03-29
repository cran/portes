\name{InvertQ}
\alias{InvertQ}

\title{Check Stationary and Invertibility of ARMA or VARMA Models}

\description{Utility function checks whether \code{ARMA} or \code{VARMA} model 
satisfies the stationary or/and the invertibility conditions.}

\usage{
InvertQ(coef)
}

\arguments{
  \item{coef}{a numeric, matrix, or array. }
}

\details{
It should be noted that, the \code{AR}(\eqn{p}) or \code{VAR}(\eqn{p}) model can always be expressed as a \eqn{kp}-dimensional 
\code{AR}(\eqn{1}) or \code{VAR}(\eqn{1}), and the \code{MA}(\eqn{q}) or \code{VMA}(\eqn{q}) model can 
always be expressed as a \eqn{kq}-dimensional \code{MA}(\eqn{1}) or \code{VMA}(\eqn{1}). 
For this reason, we can use this fact when we need to find the explicit solutions of \code{AR}(\eqn{p}) or 
\code{VAR}(\eqn{p}) models or \code{MA}(\eqn{q}) or \code{VMA}(\eqn{q}) models as the \code{AR}(\eqn{1}) or
\code{VAR}(\eqn{1}) or the \code{MA}(\eqn{1}) or \code{VMA}(\eqn{1}) models can be characterized with simple intuitive formulas.
}

\value{
A warning message only if the model is not stationary or/and not invertible.}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Lutkepohl, H. (2005). "New introduction to multiple time series analysis". 
Springer-Verlag, New York.

Reinsel, G. C. (1997). "Elements of Multivariate Time Series Analysis". 
Springer-Verlag, 2nd edition.
}

\seealso{
\code{\link{varima.sim}}, \code{\link{vma.sim}}, \code{\link{ImpulseVMA}}
}

\examples{
##############################################################
### Check Stationary
phi <- array(c(0.5,0.4,0.1,0.5,0,0.3,0,0),dim=c(2,2,2))
InvertQ(phi)
### Check Invertibility
theta <- array(c(0.5,0.4,0.1,0.5,0,0.3,0,0),dim=c(2,2,2))
InvertQ(theta)
}

\keyword{ ts}

