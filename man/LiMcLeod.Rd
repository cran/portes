\name{LiMcLeod}
\alias{LiMcLeod}

\title{The Modified Multivariate Portmanteau Test, Li-McLeod (1981)}

\description{
The modified multivariate portmanteau test suggested by Li and McLeod (1981).
}

\usage{
LiMcLeod(obj,lags=seq(5,30,5),order=0,SquaredQ=FALSE)
}

\arguments{
  \item{obj}{a univariate or multivariate series with class \code{"numeric"}, \code{"matrix"}, 
             \code{"ts"}, or \code{("mts" "ts")}.
             It can be also an object of fitted time-series model with class \code{"ar"}, 
              \code{"arima0"}, \code{"Arima"}, \code{"varest"}, \code{"FitAR"}, or \code{"FitFGN"}. }
  \item{lags}{vector of lag auto-cross correlation coefficients used for \code{LiMcLeod} test.}
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

\value{
The multivariate test statistic suggested by Li and McLeod (1981) and its corresponding p-values 
for different lags based on the asymptotic chi-square distribution with \code{k^2(lags-order)} degrees of freedom.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Li, W. K. and McLeod, A. I. (1981). "Distribution of The Residual Autocorrelations in 
Multivariate ARMA Time Series Models". Journal of The Royal Statistical Society, Series B, 43, 231-239.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link[stats]{Box.test}}, \code{\link{BoxPierce}}, \code{\link{LjungBox}}, 
\code{\link{Hosking}}, \code{\link{gvtest}}, \code{\link{portest}}, \code{\link{GetResiduals}}
}

\examples{
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
LiMcLeod(fit,lags,order = 2)        ## True
LiMcLeod(fit,lags)                  ## True
## Apply the test statistic on the residuals
res <- ts((fit$resid)[-(1:2), ])
LiMcLeod(res,lags,order = 2)        ## True
LiMcLeod(res,lags)                  ## Not True
}

\keyword{htest}
\keyword{ ts}

