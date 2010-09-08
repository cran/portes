\name{rstable}
\alias{rstable}

\title{Generate Data From Stable Distributions}

\description{
Generate data from stable distribution with infinite variance.
}

\usage{
    rstable(n, ALPHA, BETA, GAMMA = NULL, DELTA = NULL)
}

\arguments{
  \item{n}{length of the series. }
  \item{ALPHA}{index parameters, each in the range \code{(0,2]}.}
  \item{BETA}{skewness parameters, each in the range \code{[-1, 1]}.}
  \item{GAMMA}{scale parameters.}
  \item{DELTA}{location parameters.}
}

\details{
\code{ALPHA, BETA, GAMMA,} and \code{DELTA} should have the same length. 
This length, \eqn{k}, represents the number of the variables that we need to generate. 
}

\value{
A vector of dimension \eqn{n\times k} from independent stable distributions.
}

\author{ Esam Mahdi & A.I. McLeod.}

\references{
Chambers, J.M., Mallows, C.L., and Stuck, B.W. (1976). "A Method for Simulating Stable Random Variables".
Journal of American Statistical Association, 71, 340-344.
}

\seealso{
There is also a function \code{rstable} in the \code{fBasics} package for the univariate case only. 
See also \code{\link{FitStable}}, \code{\link{simvarma}}
}

\examples{
## Generate Univariate Data
n <- 500
ALPHA <- 1.75
BETA <- 0
GAMMA <- 1.5
DELTA <- 0
rstable(n, ALPHA, BETA, GAMMA, DELTA)                
## Generate Bivariate Data  
n <- 500
ALPHA <- c(1.3,1.5)
BETA <- c(0.3,-0.6)
rstable(n, ALPHA, BETA)  
## Generate Multivariate Data with 
n <- 500
ALPHA <- c(1.3,1.5,1.7)
BETA <- c(0.3,-0.6,0)
GAMMA <- c(3,1,6)
rstable(n, ALPHA, BETA,GAMMA)  
}

\keyword{distribution}







