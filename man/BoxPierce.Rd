\name{BoxPierce}
\alias{BoxPierce}

\title{The Univariate-Multivariate Box and Pierce Portmanteau Test}

\description{
The univariate or multivariate Box-Pierce (1970) portmanteau test.
}

\usage{
BoxPierce(obj,lags=seq(5,30,5),order=0,SquaredQ=FALSE)
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

\value{The Box and Pierce univariate or multivariate test statistic with the associated p-values 
for different lags based on asymptotic chi-square distribution.
}

\author{ Esam Mahdi & A.I. McLeod.}

\references{
Box, G.E.P. and Pierce, D.A. (1970). "Distribution of Residual Autocorrelation in Autoregressive-Integrated Moving Average
Time Series Models". Journal of American Statistical Association, 65, 1509-1526.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link[stats]{Box.test}}, \code{\link{LjungBox}},
 \code{\link{Hosking}}, \code{\link{LiMcLeod}}, \code{\link{GVStat}}, \code{\link{portes}}, \code{\link{Get.Resid}}
}

\examples{
x <- rnorm(100)
BoxPierce(x)
##############################################################
## Measurements of the annual flow of the river Nile at Aswan 
## from the years 1871 to 1970:
##############################################################
fit <- arima(Nile, c(2, 0, 0))
lags <- c(5, 10, 20)
## Apply the test statistic on the fitted model 
BoxPierce(fit, lags, order = 2)  ## Univariate Box-Pierce test
## Apply the test statistic on the residuals
res <- resid(fit)
BoxPierce(res, lags, order = 2)  ## Univariate Box-Pierce test
##############################################################
## Quarterly, west German investment, income, and consumption 
## from first quarter of 1960 to fourth quarter of 1982: 
##############################################################
data(WestGerman)
DiffData <- matrix(numeric(3 * 91), ncol = 3)
  for (i in 1:3) 
    DiffData[, i] <- diff(log(WestGerman[, i]), lag = 1)
fit <- ar.ols(DiffData, intercept = TRUE, order.max = 2)
lags <- c(5,10)
## Apply the test statistic on the fitted model 
BoxPierce(fit,lags,order = 2)  ## Multivariate Box-Pierce test
## Apply the test statistic on the residuals
res <- ts((fit$resid)[-(1:2), ])
BoxPierce(res,lags,order = 2)  ## Multivariate Box-Pierce test
##############################################################
## Monthly log stock returns of Intel corporation data
## Test for ARCH Effects 
##############################################################
monthintel <- as.ts(monthintel)
BoxPierce(monthintel,SquaredQ=FALSE) ## Test for usual residuals 
BoxPierce(monthintel,SquaredQ=TRUE)  ## Test for ARCH effects
}


\keyword{htest}
\keyword{ ts}

