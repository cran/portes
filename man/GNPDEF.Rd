\name{GNPDEF}
\alias{GNPDEF}
\docType{data}

\title{ 
GNP Deflator for U.S. Inflation Data from January 01, 1947 to April 01, 2010.}
\description{
GNP deflator for U.S. inflation data from 1947-01-01 to 2010-04-01.}

\usage{data(GNPDEF)}
\format{
  A data frame with 254 observations on the following 2 variables.
  \describe{
    \item{\code{time}}{time}
    \item{\code{GNPDEF}}{a numeric vector denotes the GNP deflator}
  }
}

\source{
\url{http://research.stlouisfed.org}
}

\references{
Bollerslev, T. (1986). "Generalized autoregressive conditional heteroskedasticity".
Journal of Econometrics, 31(3), 307-327.
}

\examples{
plot(ts(GNPDEF[,2]))
}

\keyword{Datasets}

