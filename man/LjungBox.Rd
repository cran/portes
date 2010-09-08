\name{LjungBox}
\alias{LjungBox}

\title{Ljung and Box Portmanteau Test}

\description{
The Ljung-Box (1978) modified portmanteau test.
}

\usage{
LjungBox(obj,lags=seq(5,30,5),order=0,SquaredQ=FALSE)
}

\arguments{
  \item{obj}{a fitted time-series model with any of the classes \code{"ar"}, \code{"Arima"}, \code{"arima0"},
   \code{"FitAR"}, or \code{"FitFGN"} or a residuals with class 
   \code{"ts"}, \code{"numeric"}, \code{"matrix"}, or \code{("mts" "ts")}. }
  \item{lags}{vector of lag values.}
  \item{order}{\eqn{p+q}, where (\eqn{p},\eqn{q}) is the \code{ARMA} or the \code{VARMA} order.}
  \item{SquaredQ}{\code{TRUE}, apply the test to the squared values. 
    This is needed to check the Autoregressive Conditional Heteroscedastic, \code{ARCH}, effects. 
    When \code{SquaredQ = FALSE}, apply for the usual residuals. }
}

\value{The Ljung and Box test statistic with the associated p-values 
for different lags based on asymptotic chi-square distribution.
}

\author{ Esam Mahdi & A.I. McLeod.}

\references{
Ljung, G.M. and Box, G.E.P (1978). "On a Measure of Lack of Fit in Time Series Models". Biometrika, 65, 297-303.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link[stats]{Box.test}}, \code{\link{BoxPierce}}, \code{\link{Hosking}}, \code{\link{LiMcLeod}},
\code{\link{GVStat}}, \code{\link{portes}}, \code{\link{Get.Resid}}
}

\examples{
x <- rnorm(100)
LjungBox(x)
##############################################################
## Measurements of the annual flow of the river Nile at Aswan 
## from the years 1871 to 1970:
##############################################################
fit <- arima(Nile, c(2, 0, 0))
lags <- c(5, 10, 20)
## Apply the test statistic on the fitted model 
LjungBox(fit, lags, order = 2)  ## Univariate LjungBox test
## Apply the test statistic on the residuals
res <- resid(fit)
LjungBox(res, lags, order = 2)  ## Univariate LjungBox test  
##############################################################
## Monthly log stock returns of Intel corporation data
## Test for ARCH Effects 
##############################################################
monthintel <- as.ts(monthintel)
LjungBox(monthintel,SquaredQ=FALSE) ## Test for usual residuals   
LjungBox(monthintel,SquaredQ=TRUE)  ## Test for ARCH effects
}

\keyword{htest}
\keyword{ ts}

