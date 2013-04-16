\name{LjungBox}
\alias{LjungBox}

\title{Ljung and Box Portmanteau Test}

\description{
The Ljung-Box (1978) modified portmanteau test.
In the multivariate time series, this test statistic is asymptotically
equal to \code{\link{Hosking}}. 
}

\usage{
LjungBox(obj,lags=seq(5,30,5),order=0,SquaredQ=FALSE)
}

\arguments{
  \item{obj}{a univariate or multivariate series with class \code{"numeric"}, \code{"matrix"}, 
             \code{"ts"}, or \code{("mts" "ts")}.
             It can be also an object of fitted time-series model with class \code{"ar"}, 
              \code{"arima0"}, \code{"Arima"}, \code{"varest"}, \code{"FitAR"}, \code{"FitFGN"},
              \code{"garch"}, or \code{"fGARCH"}. 
              \code{obj} may also an object with class \code{"list"} (see details and following example).}
  \item{lags}{vector of lag auto-cross correlation coefficients used for \code{LjungBox} test.}
  \item{order}{needed for degrees of freedom of asymptotic chi-square distribution. 
             If \code{obj} is an object with class \code{"ar"}, 
             \code{"arima0"}, \code{"Arima"}, \code{"varest"}, \code{"FitAR"}, 
             \code{"FitFGN"}, \code{"garch"}, \code{"fGARCH"}, or \code{"list"} then no need to enter 
             the value of \code{order} as it will be automatically determined.
             In general \code{order} equals to the number of estimated parameters in the fitted model.}                     
  \item{SquaredQ}{if \code{TRUE} then apply the test on the squared values. 
             This checks for Autoregressive Conditional Heteroscedastic, 
             \code{ARCH}, effects. 
             When \code{SquaredQ = FALSE}, then apply the test on the usual residuals. }
}

\details{
However the portmanteau test statistic can be applied directly on the output objects from 
the built in \code{R} functions
\code{ar()}, \code{FitAR()}, \code{arima()}, \code{arim0()}, \code{Arima()},
\code{auto.arima()}, \code{VAR()}, \code{garch()}, \code{garchFit()}, \code{FitFGN()}, etc, 
it works with output objects from any fitted model. 
In this case, users should write their own function to fit any model they want. 
The object \code{obj} represents the output of this function. 
This output must be a list with at least two outcomes: 
the fitted residual and the order of the fitted model (\code{list(res = ..., order = ...)}). 
See the following example with the function \code{FitModel()}.
}

\value{The Ljung and Box test statistic with the associated p-values 
for different lags based on the asymptotic chi-square distribution with \code{k^2(lags-order)} degrees of freedom.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Ljung, G.M. and Box, G.E.P (1978). "On a Measure of Lack of Fit in Time Series Models". 
Biometrika, 65, 297-303.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link[stats]{Box.test}}, \code{\link{BoxPierce}}, \code{\link{Hosking}}, \code{\link{LiMcLeod}},
\code{\link{gvtest}}, \code{\link{portest}}, \code{\link{GetResiduals}}, \code{\link[TSA]{tar}}
}

\examples{
x <- rnorm(100)
LjungBox(x)         
##############################################################
## Measurements of the annual flow of the river Nile at Aswan 
## from the years 1871 to 1970:
##############################################################
fit <- arima(Nile, c(1, 0, 1))
lags <- c(5, 10, 20)
## Apply the univariate test statistic on the fitted model 
LjungBox(fit, lags)               ## Correct 
LjungBox(fit, lags, order = 2)    ## Correct 
## Apply the test statistic on the residuals and set order = 2 
res <- resid(fit)
LjungBox(res, lags)               ## Wrong  
LjungBox(res, lags, order = 2)    ## Correct 
#to save time the other examples are not run
\dontrun{
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
LjungBox(fit,lags)                ## Correct 
## Apply the test statistic on the residuals where order = 2
res <- ts((fit$resid)[-(1:2), ])
LjungBox(res,lags)                ## Wrong 
LjungBox(res,lags,order = 2)      ## Correct 
##############################################################
## Monthly log stock returns of Intel corporation data
## Test for ARCH Effects 
##############################################################
monthintel <- as.ts(monthintel)
LjungBox(monthintel)                ## Usual test 
LjungBox(monthintel,SquaredQ=TRUE)  ## Test for ARCH effects
##############################################################
## Write a function to fit a model 
## Apply portmanteau test on fitted obj with class "list"
##############################################################
## Example 1
library("FitAR")
FitModel <- function(data){
    fit <- FitAR(z=data,p=2)
    p <- length(fit$phiHat)
    order <- p
    res <- fit$res 
 list(res=res,order=order)
}
Fit <- FitModel(Nile)
LjungBox(Fit)
##
## Example 2
library("TSA")
FitModel <- function(data){
    fit <- TSA::tar(y=log(data),p1=4,p2=4,d=3,a=0.1,b=0.9,print=FALSE)
    res <- ts(fit$std.res)
    p1 <- fit$p1
    p2 <- fit$p2
    order <- max(p1, p2)
    parSpec <- list(res=res,order=order)
  parSpec
}
data(prey.eq)
Fit <- FitModel(prey.eq) 
LjungBox(Fit)
}
}

\keyword{htest}
\keyword{ ts}

