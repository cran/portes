\name{EconomicUK}
\alias{EconomicUK}
\docType{data}
\title{
Quarterly U.K. economic time series from 1957 Q3 to 1967 Q4
}
\description{
The data are quarterly, seasonally unadjusted in 1958 prices, covering the period
1957/3-1967/4 (with 7 series each with 42 observations), as published in Economic Trends, with information 
about consumers' expenditure on goods and services, Investment, inventory investment, 
imports of goods and services, gross domestic product, and personal disposable income.
Prothero and Wallis (1976) fitted several models to each series and compared
their performance with a multivariate model.
}
\usage{data("EconomicUK")}
\format{
  A data frame with 42 observations on the following 8 variables.
  \describe{
    \item{\code{Cd}}{consumers' expenditure on durable goods}
    \item{\code{Cn}}{consumers' expenditure on all other goods and services}
    \item{\code{I}}{investment (gross domestic fixed capital formation)}
    \item{\code{Iv}}{inventory investment (value of physical increase in stocks and work in progress)}
    \item{\code{M}}{imports of goods and services}
    \item{\code{Y}}{gross domestic product}
    \item{\code{Yd}}{personal disposable income}
    \item{\code{year}}{year with attributed number associated to quarterly period}
  }
}

\source{
The data are quarterly, seasonally unadjusted in 1958 prices, covering the period
1957/3-1967/4 (42 observations), as published in Economic Trends.
}

\references{
David L. Prothero and Kenneth F. Wallis (1976). "Modelling macroeconomic time series (with discussion)", Journal of the Royal Statistical Society, A, Vol.139, Part 4, pp.468-500.
}

\keyword{datasets}
