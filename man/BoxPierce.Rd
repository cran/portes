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
  \item{obj}{a univariate or multivariate series with class \code{"numeric"}, \code{"matrix"}, 
             \code{"ts"}, or \code{("mts" "ts")}.
             It can be also an object of fitted time-series model with class \code{"ar"}, 
              \code{"arima0"}, \code{"Arima"}, \code{"varest"}, \code{"FitAR"}, or \code{"FitFGN"}. }
  \item{lags}{vector of lag auto-cross correlation coefficients used for \code{BoxPierce} test.}
  \item{order}{this is used for degrees of freedom of asymptotic chi-square distribution. 
             If \code{obj} is a fitted time-series model with class \code{"ar"}, 
             \code{"arima0"}, \code{"Arima"}, \code{"varest"}, \code{"FitAR"}, or \code{"FitFGN"} then no need to enter 
             the value of \code{order} as it will be automatically determined.
             In general \code{order = p + q}, where \code{p} and \code{q} are the orders of the autoregressive (or vector autoregressive) 
             and moving average (or vector moving average) models respectively. 
             \code{order = 0} is used for testing random series, fractional gaussian noise, 
             or generalized autoregressive conditional heteroscedasticity.}           
  \item{SquaredQ}{if \code{TRUE} then apply the test on the squared values. 
             This checks for Autoregressive Conditional Heteroscedastic, \code{ARCH}, effects. 
             When \code{SquaredQ = FALSE}, then apply the test on the usual residuals. }
}

\value{The Box and Pierce univariate or multivariate test statistic with the associated p-values 
for different lags based on the asymptotic chi-square distribution.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Box, G.E.P. and Pierce, D.A. (1970). "Distribution of Residual Autocorrelation in Autoregressive-Integrated Moving Average
Time Series Models". Journal of American Statistical Association, 65, 1509-1526.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link[stats]{Box.test}}, \code{\link{LjungBox}},
 \code{\link{Hosking}}, \code{\link{LiMcLeod}}, 
\code{\link{gvtest}}, \code{\link{portest}}, \code{\link{GetResiduals}}
}

\examples{
x <- rnorm(100)
BoxPierce(x)         
##############################################################
## Measurements of the annual flow of the river Nile at Aswan 
## from the years 1871 to 1970:
##############################################################
fit <- arima(Nile, c(1, 0, 1))
lags <- c(5, 10, 20)
## Apply the univariate test statistic on the fitted model 
BoxPierce(fit, lags)               ## True 
BoxPierce(fit, lags, order = 2)    ## True 
## Apply the test statistic on the residuals and set order = 2 
res <- resid(fit)
BoxPierce(res, lags)               ## Not True  
BoxPierce(res, lags, order = 2)    ## True 
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
BoxPierce(fit,lags)                ## True 
## Apply the test statistic on the residuals where order = 2
res <- ts((fit$resid)[-(1:2), ])
BoxPierce(res,lags)                ## Not True 
BoxPierce(res,lags,order = 2)      ## True 
##############################################################
## Monthly log stock returns of Intel corporation data
## Test for ARCH Effects 
##############################################################
monthintel <- as.ts(monthintel)
BoxPierce(monthintel)                ## Usual test 
BoxPierce(monthintel,SquaredQ=TRUE)  ## Test for ARCH effects
}

\keyword{htest}
\keyword{ ts}

