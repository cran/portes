\name{FitStable}
\alias{FitStable}

\title{ Fit Parameters to Stable Distributions, McCulloch (1986) }

\description{
The quantile method of McCullogh (1986).
}

\usage{
FitStable(x)
}

\arguments{
  \item{x}{ univariate or independent multivariate data. }
}

\details{
The quantile estimation method of McCulloch (1986) is used for each variable in x.
It is highly reliable, fast and reasonably efficient especially bearing in
mind that in most applications there is a lot of data.  
}

\value{
k rows represents the number of variables in the vector x, and 4 columns with named components 
alpha, beta, scale, and location respectively. 
}

\author{ A.I. McLeod and Jen-Wen Lin}

\references{
Lin, Jen-Wen and McLeod (2008). 
"Portmanteau Tests for ARMA Models with Infinite Variance."
Journal of Time Series Analysis, 29, 600-617.
 
McCulloch, J. H. (1986). "Simple Consistent Estimator of Stable Distribution Parameters".
Commun. Statist.--Simula., 15(4), 1109-1136.}

\seealso{There is also a function \code{stableFit} in the \code{fBasics} package for fitting
stable distributions for univariate data but it does not work very well. This is the reason
for \code{FitStable}. See also \code{\link{rstable}}, \code{\link{simvarma}}, 
}

\examples{
## Univariate
x <- rstable(800, 1.7, 0, 1, 0)
FitStable(x)
## Multivariate
ALPHA <- c(1.3,1.6)
BETA <- c(0,0.2)
GAMMA <-c(1,1)
DELTA <-c(0,0.2)
x <- rstable(500, ALPHA, BETA, GAMMA, DELTA)
FitStable(x)
}

\keyword{distribution}
