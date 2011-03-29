\name{CRSP}
\alias{CRSP}
\docType{data}

\title{ Monthly simple returns of the CRSP value-weighted index, 1926 to 1997}

\description{
This data has been discussed by Tsay (2002, Ch.2, p.38 and 39) and
Lin and McLeod (2008). 
There are 864 data values.
}

\usage{data(CRSP)}

\references{
Lin, J.-W. and McLeod, A.I. (2008). "Portmanteau Tests for ARMA Models with Infinite Variance".  
Journal of Time Series Analysis, 29, 600-617.

Tsay, R. S. (2002). Analysis of Financial Time Series. Wiley, New York.
}

\examples{
acf(CRSP)
fitstable(CRSP)
}

\keyword{datasets}
