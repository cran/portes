\name{GVStat}
\alias{GVStat}

\title{Generalized Variance Portmanteau Test}

\description{
New generalized variance portmanteau test, based on the determinant of the Hosking's autocorrelation block 
Toeplitz matrix with order \eqn{m+1} given in the function \code{\link{blockToeplitz}}, 
where \eqn{m} represents the order of the block matrix. 
}

\usage{
GVStat(obj,lags=seq(5,30,5),order=0,SquaredQ=FALSE)
}

\arguments{
  \item{obj}{a fitted time-series model with any of the classes \code{"ar"}, \code{"Arima"}, \code{"arima0"},
   \code{"FitAR"}, or \code{"FitFGN"} or a residuals with class 
   \code{"ts"}, \code{"numeric"}, \code{"matrix"}, or \code{("mts" "ts")}. }
  \item{lags}{vector of lag values. }
  \item{order}{\eqn{p+q}, where (\eqn{p},\eqn{q}) is the \code{ARMA} or the \code{VARMA} order.}
  \item{SquaredQ}{\code{TRUE}, apply the test to the squared values. 
    This is needed to check the Autoregressive Conditional Heteroscedastic, \code{ARCH}, effects. 
    When \code{SquaredQ = FALSE}, apply for the usual residuals. }
}

\value{
The test statistic and its associated p-values for different lags based on asymptotic 
chi-square distribution with \code{k^2(lags-order)} degrees of freedom.
}

\author{ Esam Mahdi & A.I. McLeod.}

\references{
Pena, D. and Rodriguez, J. (2002). "A Powerful Portmanteau Test of Lack of Test for Time Series".
Journal of American Statistical Association, 97, 601-610.

Mahdi, E. and McLeod, A.I. (2010). "Improved multivariate portmanteau diagnostic test". Submitted.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link[stats]{Box.test}}, \code{\link{BoxPierce}}, \code{\link{LjungBox}}, \code{\link{Hosking}}, 
\code{\link{LiMcLeod}}, \code{\link{blockToeplitz}}, \code{\link{portes}}, \code{\link{Get.Resid}}
}


\examples{
##############################################################
## Measurements of the annual flow of the river Nile at Aswan 
## from the years 1871 to 1970:
##############################################################
fit <- arima(Nile, c(2, 0, 0))
lags <- c(5, 10, 20)
## Apply the test statistic on the fitted model 
GVStat(fit, lags, order = 2)  ## Univariate GVStat test
## Apply the test statistic on the residuals
res <- resid(fit)
GVStat(res, lags, order = 2)  ## Univariate GVStat test
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
GVStat(fit,lags,order = 2)  ## Multivariate GVStat test
## Apply the test statistic on the residuals
res <- ts((fit$resid)[-(1:2), ])
GVStat(res,lags,order = 2)  ## Multivariate GVStat test
##############################################################
## Monthly log stock returns of Intel corporation data
## Test for ARCH Effects 
##############################################################
monthintel <- as.ts(monthintel)
GVStat(monthintel,SquaredQ=FALSE) ## Test for usual residuals   
GVStat(monthintel,SquaredQ=TRUE)  ## Test for ARCH effects
}

\keyword{htest}
\keyword{ ts}

