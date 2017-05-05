\name{ToeplitzBlock}
\alias{ToeplitzBlock}

\title{Toeplitz Block Matrix of Hosking (1980) Auto and Cross Correlation Matrices}

\description{Block Toeplitz matrix of order \eqn{m+1} with \eqn{k\times k} auto-cross correlation matrices. 
The Hosking (1980) definition of the correlation matrix is used.
This is needed for the function \code{\link{MahdiMcLeod}}. 
}

\usage{
ToeplitzBlock(res,lag.max,season=1)
}

\arguments{
  \item{res}{residuals, numeric or matrix. }
  \item{lag.max}{an integer number = \eqn{m} is used to determined the order of the block matrix. }
  \item{season}{seasonal periodicity given from \code{\link{MahdiMcLeod}}. Default is 1 for non seasonality cases.}
}

\value{ 
A block Toeplitz matrix of auto and cross correlation matrices using Hosking (1980) definition 
from \code{lag} = \eqn{0} to \code{lag} = \eqn{m}.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Hosking, J. R. M. (1980). "The Multivariate Portmanteau Statistic". 
Journal of American Statistical Association, 75, 602-608.

Lin, J.-W. and McLeod, A.I. (2006). "Improved Generalized Variance Portmanteau Test". 
Computational Statistics and Data Analysis, 51, 1731-1738.

Mahdi, E. and McLeod, A.I. (2011, accepted). "Improved Multivariate Portmanteau Test". 
Journal of Time Series Analysis. (JTSA - 3192).
}

\seealso{
\code{\link[stats]{acf}}, \code{\link{MahdiMcLeod}}, \code{\link[stats]{toeplitz}}
}

\examples{
x <- rnorm(100)   
ToeplitzBlock(x,lag.max=4)          ## Univariate Series
#
y <- cbind(rnorm(100),rnorm(100)) 
ToeplitzBlock(y,lag.max=4)          ## Multivariate Series
#
ToeplitzBlock(y,lag.max=4,season=4) ## Multivariate Series
}

\keyword{Array}

