\name{monthibmspln}
\alias{monthibmspln}
\docType{data}

\title{The Monthly Log Returns of IBM Stock and the S&P 500 Index}

\description{
The monthly log returns of IBM stock and the S&P 500 index from January 1926 to December 1999. 
This data has been discussed by Tsay (2005, p.356).
}

\usage{data(monthibmspln)}
\format{
  A data frame with 888 observations on the following 2 variables.
  \describe{
    \item{\code{IBM}}{a numeric vector}
    \item{\code{SP500}}{a numeric vector}
  }
}

\source{
\url{http://faculty.chicagobooth.edu/ruey.tsay/teaching/fts2/}
}

\references{
Tsay, R. S. (2005). "Analysis of Financial Time Series". Wiley, New York, 2nd edition.
}

\examples{
acf(monthibmspln)
fitstable(monthibmspln)
}

\keyword{datasets}

