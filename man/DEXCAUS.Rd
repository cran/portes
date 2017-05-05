\name{DEXCAUS}
\alias{DEXCAUS}
\docType{data}

\title{ Canada/US Foreign Exchanges Rates, Daily, Jan. 4, 1971 to Sept. 5, 1996. }

\description{
There are 2513 data values.
}

\usage{data(DEXCAUS)}

\details{
Title:               Canada / U.S. Foreign Exchange Rate
Series ID:           DEXCAUS
Source:              Board of Governors of the Federal Reserve System
Release:             H.10 Foreign Exchange Rates
Seasonal Adjustment: Not Applicable
Frequency:           Daily
Units:               Canadian Dollars to One U.S. Dollar
Date Range:          1971-01-04 to 2006-09-05
Last Updated:        2006-09-06 10:42 AM CT
Notes:               Noon buying rates in New York City for cable transfers payable in foreign currencies.
}

\source{
\url{http://research.stlouisfed.org/fred2/series/DEXCAUS}
}

\examples{
acf(DEXCAUS)
fitstable(DEXCAUS)
}

\keyword{Datasets}

