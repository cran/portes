\name{IbmSp500}
\alias{IbmSp500}
\docType{data}
\title{Monthly Returns of IBM and S&P 500 Index
}
\description{
The monthly returns of IBM stock and the S&P 500 index from January 1926 to December 2008. 
This data has been discussed by Tsay (2010, Chapter 8).
}
\usage{data(IbmSp500)}
\format{
  A data frame with 996 observations on the following 3 variables.
  \describe{
    \item{\code{date}}{a numeric vector}
    \item{\code{ibm}}{a numeric vector}
    \item{\code{sp}}{a numeric vector}
  }
}
\source{
http://faculty.chicagobooth.edu/ruey.tsay/teaching/fts3/}
\references{
Tsay, R. S. (2010). "Analysis of Financial Time Series". Wiley, New York, 3rd edition.
}
\examples{
data(IbmSp500)
plot(IbmSp500) 
acf(IbmSp500)
fitstable(IbmSp500)

}
\keyword{Datasets}
