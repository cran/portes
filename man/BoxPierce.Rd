\name{BoxPierce}
\alias{BoxPierce}

\title{The Univariate-Multivariate Box and Pierce Portmanteau Test}

\description{
The univariate or multivariate Box-Pierce (1970) portmanteau test.
}

\usage{
BoxPierce(obj,lags=seq(5,30,5),order=0,season=1,squared.residuals=FALSE)
}

\arguments{
  \item{obj}{a univariate or multivariate series with class \code{"numeric"}, \code{"matrix"}, 
             \code{"ts"}, or \code{("mts" "ts")}.
             It can be also an object of fitted time-series model with class \code{"ar"}, 
              \code{"arima0"}, \code{"Arima"}, \code{("ARIMA forecast ARIMA Arima")}, \code{"lm"},
              \code{("glm" "lm")}, or \code{"varest"}. 
              \code{obj} may also an object with class \code{"list"} (see details and following examples).}
  \item{lags}{vector of lag auto-cross correlation coefficients used for \code{Hosking} test.}
  \item{order}{Default is zero for testing the randomness of a given sequence with 
              class \code{"numeric"}, \code{"matrix"}, \code{"ts"}, or \code{("mts" "ts")}. 
              In general \code{order} equals to the number of estimated parameters in the fitted model.
             If \code{obj} is an object with class \code{"ar"}, \code{"arima0"}, \code{"Arima"}, 
             \code{"varest"}, \code{("ARIMA forecast ARIMA Arima")}, or \code{"list"} then no need to enter 
             the value of \code{order} as it will be automatically determined.
             For \code{obj} with other classes, the \code{order} is needed for degrees 
             of freedom of asymptotic chi-square distribution.}           
  \item{season}{seasonal periodicity for testing seasonality. Default is 1 for testing the non seasonality cases.}
  \item{squared.residuals}{if \code{TRUE} then apply the test on the squared values. 
             This checks for Autoregressive Conditional Heteroscedastic, 
             \code{ARCH}, effects. 
             When \code{squared.residuals = FALSE}, then apply the test on the usual residuals. }
}

\details{
However the portmanteau test statistic can be applied directly on the output objects from 
the built in \code{R} functions \code{ar()}, \code{ar.ols()}, \code{ar.burg()}, 
\code{ar.yw()}, \code{ar.mle()}, \code{arima()}, \code{arim0()}, \code{Arima()},
\code{auto.arima()}, \code{lm()}, \code{glm()}, and \code{VAR()}, 
it works with output objects from any fitted model. 
In this case, users should write their own function to fit any model they want, where they 
may use the built in \code{R} functions \code{FitAR()}, \code{garch()}, \code{garchFit()}, 
\code{fracdiff()}, \code{tar()}, etc.
The object \code{obj} represents the output of this function. 
This output must be a list with at least two outcomes: 
the fitted residual and the order of the fitted model (\code{list(res = ..., order = ...)}). 
See the following example with the function \code{FitModel()}.

Note: In \code{stats} \code{R}, the function \code{\link[stats]{Box.test}} was built to compute
the Box and Pierce (1970) and Ljung and Box (1978) test statistics
only in the univariate case where we can not use more than one single lag value at a time.
The functions \code{\link{BoxPierce}} and \code{\link{LjungBox}} are more accurate than
\code{\link[stats]{Box.test}} function and can be used in the univariate or multivariate time series
at vector of different lag values as well as they can be applied on an output object 
from a fitted model described in the description of the function \code{\link{BoxPierce}}.
}

\value{The Box and Pierce univariate or multivariate test statistic with the associated p-values 
for different lags based on the asymptotic chi-square distribution with \code{k^2(lags-order)} degrees of freedom.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Box, G.E.P. and Pierce, D.A. (1970). "Distribution of Residual Autocorrelation in Autoregressive-Integrated Moving Average
Time Series Models". Journal of American Statistical Association, 65, 1509-1526.
}

\seealso{
\code{\link[stats]{Box.test}}, \code{\link{LjungBox}}, \code{\link{MahdiMcLeod}}, \code{\link{Hosking}}, 
\code{\link{LiMcLeod}}, \code{\link{portest}}, \code{\link{GetResiduals}}.
}

\examples{
x <- rnorm(100)
BoxPierce(x)                              ## univariate test
x <- cbind(rnorm(100),rnorm(100))
BoxPierce(x)                              ## multivariate test      
##
##
## Annual flow of the river Nile at Aswan - 1871 to 1970
fit <- arima(Nile, c(1, 0, 1))
lags <- c(5, 10, 20)
## Apply the univariate test statistic on the fitted model 
BoxPierce(fit, lags)            ## Correct (no need to specify order) 
BoxPierce(fit, lags, order = 2) ## Correct 
## Apply the test statistic on the residuals and set order = 2 
res <- resid(fit)
BoxPierce(res, lags)             ## Wrong (order is needed!)  
BoxPierce(res, lags, order = 2)  ## Correct 
##
##
## Quarterly, west German investment, income, and consumption from 1960 Q1 to 1982 Q4 
data(WestGerman)
DiffData <- matrix(numeric(3 * 91), ncol = 3)
  for (i in 1:3) 
    DiffData[, i] <- diff(log(WestGerman[, i]), lag = 1)
fit <- ar.ols(DiffData, intercept = TRUE, order.max = 2)
lags <- c(5,10)
## Apply the test statistic on the fitted model 
BoxPierce(fit,lags)                ## Correct (no need to specify order)
## Apply the test statistic on the residuals where order = 2
res <- ts((fit$resid)[-(1:2), ])
BoxPierce(res,lags)                ## Wrong (order is needed!)  
BoxPierce(res,lags,order = 2)      ## Correct 
##
##
## Monthly log stock returns of Intel corporation data: Test for ARCH Effects 
monthintel <- as.ts(monthintel)
BoxPierce(monthintel)                         ## Usual test 
BoxPierce(monthintel,squared.residuals=TRUE)  ## Test for ARCH effects
##
##
## Test for seasonality
## Accidental Deaths in the US 1973 - 1978
seasonal.arima <- arima(USAccDeaths, order = c(0,1,1), seasonal = list(order = c(0,1,1)))
BoxPierce(seasonal.arima, lags = 5, season = 12)
## Quarterly U.K. economic time series from 1957 Q3 to 1967 Q4
cd <- EconomicUK[,1]
cd.fit <- arima(cd,order=c(0,1,0),seasonal=list(order=c(0,1,1),period=4))
BoxPierce(cd.fit, lags = c(5,10), season = 4)
##
##
#### Write a function to fit a model: Apply portmanteau test on fitted obj with class "list"
## Example 1 
require("FitAR")
FitModel <- function(data){
    fit <- FitAR(z=data,p=2)
    p <- length(fit$phiHat)
    order <- p
    res <- fit$res 
 list(res=res,order=order)
}
Fit <- FitModel(Nile)
BoxPierce(Fit) 
detach(package:FitAR)
##
## Example 2
FitModel <- function(data){
    fit <- ar.ols(data, intercept = TRUE, order.max = 2)
    order <- 2
    res <- res <- ts((fit$resid)[-(1:2), ]) 
 list(res=res,order=order)
}
data(WestGerman)
DiffData <- matrix(numeric(3 * 91), ncol = 3)
  for (i in 1:3) 
    DiffData[, i] <- diff(log(WestGerman[, i]), lag = 1)
Fit <- FitModel(DiffData)
BoxPierce(Fit) 
}

\keyword{Portmanteau Test}


