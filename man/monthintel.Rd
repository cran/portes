\name{monthintel}
\alias{monthintel}
\docType{data}

\title{The Monthly Log Stock Returns of Intel Corporation from January 1973 to December 2003
}

\description{
The monthly log stock returns of Intel Corporation from January 1973 to December 2003.
This data has been discussed by Tsay (2005, p.99-102).
There are 372 data values.
}

\usage{data(monthintel)}

\source{
\url{http://faculty.chicagobooth.edu/ruey.tsay/teaching/fts2/}
}

\references{
Tsay, R. S. (2005). "Analysis of Financial Time Series". Wiley, New York, 2nd edition.
}

\examples{
acf(monthintel)
fitstable(monthintel)
}

\keyword{Datasets}

