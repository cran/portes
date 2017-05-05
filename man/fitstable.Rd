\name{fitstable}
\alias{fitstable}

\title{ Fit Parameters to Stable Distributions, McCulloch (1986) }

\description{
The quantile method of McCullogh (1986).
}

\usage{
fitstable(x)
}

\arguments{
  \item{x}{ univariate or independent multivariate variables of dimension k. }
}

\details{
The quantile estimation method of McCulloch (1986) is used for each variable in x.
It is highly reliable, fast and reasonably efficient especially bearing in
mind that in most applications there is a lot of data.  
}

\value{
matrix of k rows and 4 columns. 
k represents the number of the variables in the vector x and the columns with named components 
alpha, beta, scale, and location respectively. 
}

\author{Esam Mahdi, A.I. McLeod, and Jen-Wen Lin.}

\references{
Lin, J.-W. and McLeod A.I.(2008). 
"Portmanteau Tests for ARMA Models with Infinite Variance."
Journal of Time Series Analysis, 29, 600-617.
 
McCulloch, J. H. (1986). "Simple Consistent Estimator of Stable Distribution Parameters".
Commun. Statist.--Simula., 15(4), 1109-1136.}

\seealso{There is also a function \code{stableFit()} in the \code{fBasics} package for fitting
stable distributions for univariate data. 
See also \code{\link{rStable}}, \code{\link{varima.sim}}, 
}

\examples{
## Univariate
x <- rStable(800, 1.7, 0, 1, 0)
fitstable(x)
## Multivariate
ALPHA <- c(1.3,1.6)
BETA <- c(0,0.2)
GAMMA <-c(1,1)
DELTA <-c(0,0.2)
x <- rStable(500, ALPHA, BETA, GAMMA, DELTA)
fitstable(x)
}

\keyword{Distribution}
