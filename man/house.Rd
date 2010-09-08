\name{house}
\alias{house}
\docType{data}

\title{
The Monthly House Sales and House Starts.
}

\description{
The monthly house sales (in thousands of units) and house starts (in thousands of units) for the time period January 1965 and December 1975.
This data has been discussed by Wei (1990, p.330) and
McQuarrie and Tsai (1998, p.409). 
}

\usage{data(house)}

\format{
  A data frame with 132 observations on the following 2 variables.
  \describe{
    \item{\code{sale}}{house sales (in thousands)}
    \item{\code{start}}{house starts (in thousands)}
  }
}

\source{
Wei, W. W. S., 1990. "Time Series Analysis". Addison-Wesley, New York, 1st edition. 
}

\references{
Wei, W. W. S., 1990. "Time Series Analysis". Addison-Wesley, New York, 1st edition.

McQuarrie, A. D., Tsai, C.-L., 1998. "Regression and time series model selection". World Scientific.
}

\examples{
acf(CRSP)
FitStable(CRSP)
}

\keyword{datasets}

