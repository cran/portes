\name{ToeplitzBlock}
\alias{ToeplitzBlock}

\title{Toeplitz Block Matrix of Hosking (1980) Auto and Cross Correlation Matrices}

\description{Block Toeplitz matrix of order \eqn{m+1} with \eqn{k\times k} auto-cross correlation matrices. 
The Hosking (1980) definition of the correlation matrix is used.
This function is needed for the function \code{\link{gvtest}}. 
}

\usage{
ToeplitzBlock(res,lag.max)
}

\arguments{
  \item{res}{residuals, numeric or matrix. }
  \item{lag.max}{an integer number = \eqn{m} is used to determined the order of the block matrix. }
}

\value{ 
A block Toeplitz matrix of auto and cross correlation matrices based on Hosking (1980) definition 
from \code{lag} = \eqn{0} to \code{lag} = \eqn{m}.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Hosking, J. R. M. (1980). "The Multivariate Portmanteau Statistic". 
Journal of American Statistical Association, 75, 602-608.

Lin, J.-W. and McLeod, A.I. (2006). "Improved Generalized Variance Portmanteau Test". 
Computational Statistics and Data Analysis, 51, 1731-1738.

Mahdi, E. and McLeod, A.I. (2011). "Improved multivariate portmanteau diagnostic test". Submitted.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link{gvtest}}, \code{\link[stats]{toeplitz}}
}

\examples{
## Univariate Series
x <- rnorm(100)
ToeplitzBlock(x,lag.max=4)
## Multivariate Series
x <- cbind(rnorm(100),rnorm(100))
ToeplitzBlock(x,lag.max=4)
}

\keyword{array}

