\name{gvtest}
\alias{gvtest}

\title{Generalized Variance Portmanteau Test}

\description{
New generalized variance portmanteau test based on the determinant of the Hosking's autocorrelation block 
Toeplitz matrix with order \eqn{m+1} given in the function \code{\link{ToeplitzBlock}}, 
where \eqn{m} represents the order of the block matrix. 
}

\usage{
gvtest(obj,lags=seq(5,30,5),order=0,SquaredQ=FALSE)
}

\arguments{
  \item{obj}{a univariate or multivariate series with class \code{"numeric"}, \code{"matrix"}, 
             \code{"ts"}, or \code{("mts" "ts")}.
             It can be also an object of fitted time-series model with class \code{"ar"}, 
              \code{"arima0"}, \code{"Arima"}, \code{"varest"}, \code{"FitAR"}, or \code{"FitFGN"}. }
  \item{lags}{vector of lag auto-cross correlation coefficients used for \code{gvtest} test.}
  \item{order}{this is used for degrees of freedom of asymptotic gamma distribution. 
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
The generalized variance portmanteau test statistic and its associated p-values for different lags based on asymptotic 
gamma distribution as suggested by Pena and Rodriguez (2002, 2006) in the univariate time series and generalized to
the multivariate time series by Mahdi and McLeod (2011).
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Mahdi, E. and McLeod, A.I. (2011). "Improved multivariate portmanteau diagnostic test". Submitted.

Pena, D. and Rodriguez, J. (2002). "A Powerful Portmanteau Test of Lack of Test for Time Series".
Journal of American Statistical Association, 97, 601-610.

Pena, D. and Rodriguez, J. (2006). "The log of the determinant of the autocorrelation
matrix for testing goodness of fit in time series". Journal of Statistical Planning and
Inference, 136, 2706-2718.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link[stats]{Box.test}}, \code{\link{BoxPierce}}, 
\code{\link{LjungBox}}, \code{\link{Hosking}}, 
\code{\link{LiMcLeod}}, \code{\link{portest}}, \code{\link{ToeplitzBlock}}, \code{\link{GetResiduals}}
}

\examples{
x <- rnorm(100)
gvtest(x)            
##############################################################
## Measurements of the annual flow of the river Nile at Aswan 
## from the years 1871 to 1970:
##############################################################
fit <- arima(Nile, c(1, 0, 1))
lags <- c(5, 10, 20, 30)
## Apply the univariate test statistic on the fitted model 
gvtest(fit, lags)               ## True 
gvtest(fit, lags, order = 2)    ## True 
## Apply the test statistic on the residuals and set order = 2 
res <- resid(fit)
gvtest(res, lags)               ## Not True  
gvtest(res, lags, order = 2)    ## True 
##############################################################
## Quarterly, west German investment, income, and consumption 
## from first quarter of 1960 to fourth quarter of 1982: 
##############################################################
data(WestGerman)
DiffData <- matrix(numeric(3 * 91), ncol = 3)
  for (i in 1:3) 
    DiffData[, i] <- diff(log(WestGerman[, i]), lag = 1)
fit <- ar.ols(DiffData, intercept = TRUE, order.max = 2)
lags <- seq(5,30,5)
## Apply the test statistic on the fitted model 
gvtest(fit,lags)                ## True 
## Apply the test statistic on the residuals where order = 2
res <- ts((fit$resid)[-(1:2), ])
gvtest(res,lags)                ## Not True 
gvtest(res,lags,order = 2)      ## True 
##############################################################
## Monthly log stock returns of Intel corporation data
## Test for ARCH Effects 
##############################################################
monthintel <- as.ts(monthintel)
gvtest(monthintel)                ## Usual test
gvtest(monthintel,SquaredQ=TRUE)  ## Test for ARCH effects
}

\keyword{htest}
\keyword{ ts}

