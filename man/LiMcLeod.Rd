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
  \item{obj}{a fitted time-series model with any of the classes \code{"ar"}, \code{"Arima"}, \code{"arima0"},
   \code{"FitAR"}, or \code{"FitFGN"} or a residuals with class 
   \code{"ts"}, \code{"numeric"}, \code{"matrix"}, or \code{("mts" "ts")}. }
  \item{lags}{vector of lag values.}
  \item{order}{\eqn{p+q}, where (\eqn{p},\eqn{q}) is the \code{ARMA} or the \code{VARMA} order.}
  \item{SquaredQ}{\code{TRUE}, apply the test to the squared values. 
    This is needed to check the Autoregressive Conditional Heteroscedastic, \code{ARCH}, effects. 
    When \code{SquaredQ = FALSE}, apply for the usual residuals. }
}

\value{
The multivariate test statistic suggested by Li and McLeod (1981) and its corresponding p-values 
for different lags based on asymptotic chi-square distribution with \code{k^2(lags-order)} degrees of freedom.
}

\author{ Esam Mahdi & A.I. McLeod.}

\references{
Li, W. K. and McLeod, A. I. (1981). "Distribution of The Residual Autocorrelations in 
Multivariate ARMA Time Series Models". Journal of The Royal Statistical Society, Series B, 43, 231-239.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link[stats]{Box.test}}, \code{\link{BoxPierce}}, \code{\link{LjungBox}}, 
\code{\link{Hosking}},\code{\link{GVStat}}, \code{\link{portes}}, \code{\link{Get.Resid}}
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
LiMcLeod(fit,lags,order = 2)  
## Apply the test statistic on the residuals
res <- ts((fit$resid)[-(1:2), ])
LiMcLeod(res,lags,order = 2)  
}

\keyword{htest}
\keyword{ ts}

