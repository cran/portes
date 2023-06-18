\name{LjungBox}
\alias{LjungBox}

\title{Ljung and Box Portmanteau Test}

\description{
The Ljung-Box (1978) modified portmanteau test.
In the multivariate time series, this test statistic is asymptotically
equal to \code{\link{Hosking}}. 
}

\usage{
LjungBox(obj,lags=seq(5,30,5),fitdf=0,sqrd.res=FALSE)
}

\arguments{
  \item{obj}{a univariate or multivariate series with class \code{"numeric"}, \code{"matrix"}, 
             \code{"ts"}, or \code{("mts" "ts")}.
             It can be also an object of fitted time-series model with class \code{"ar"}, 
              \code{"arima0"}, \code{"Arima"}, \code{("ARIMA forecast ARIMA Arima")}, \code{"lm"},
              \code{("glm" "lm")}, or \code{"varest"}. 
              \code{obj} may also an object with class \code{"list"} (see details and following examples).}
  \item{lags}{vector of lag auto-cross correlation coefficients used for \code{Hosking} test.}
  \item{fitdf}{Default is zero for testing the randomness of a given sequence with 
              class \code{"numeric"}, \code{"matrix"}, \code{"ts"}, or \code{("mts" "ts")}. 
              In general \code{fitdf} equals to the number of estimated parameters in the fitted model.
             If \code{obj} is an object with class \code{"ar"}, \code{"arima0"}, \code{"Arima"}, 
             \code{"varest"}, \code{("ARIMA forecast ARIMA Arima")}, or \code{"list"} then no need to enter 
             the value of \code{fitdf} as it will be automatically determined.
             For \code{obj} with other classes, the \code{fitdf} is needed for degrees 
             of freedom of asymptotic chi-square distribution.}           
  \item{sqrd.res}{if \code{TRUE} then apply the test on the squared values. 
             This checks for Autoregressive Conditional Heteroscedastic, 
             \code{ARCH}, effects. 
             When \code{sqrd.res = FALSE}, then apply the test on the usual residuals. }
}

\details{
However the portmanteau test statistic can be applied directly on the output objects from 
the built in \code{R} functions \code{ar()}, \code{ar.ols()}, \code{ar.burg()}, 
\code{ar.yw()}, \code{ar.mle()}, \code{arima()}, \code{arim0()}, \code{Arima()},
\code{auto.arima()}, \code{lm()}, \code{glm()}, and \code{VAR()}, 
it works with output objects from any fitted model. 
In this case, users should write their own function to fit any model they want, where they 
may use the built in \code{R} functions \code{garch()}, \code{garchFit()}, 
\code{fracdiff()}, \code{tar()}, etc.
The object \code{obj} represents the output of this function. 
This output must be a list with at least two outcomes: 
the fitted residual and the fitdf of the fitted model (\code{list(res = ..., fitdf = ...)}). 
See the following example with the function \code{FitModel()}.

Note: In \code{stats} \code{R}, the function \code{\link[stats]{Box.test}} was built to compute
the Box and Pierce (1970) and Ljung and Box (1978) test statistics
only in the univariate case where we can not use more than one single lag value at a time.
The functions \code{\link{BoxPierce}} and \code{\link{LjungBox}} are more accurate than
\code{\link[stats]{Box.test}} function and can be used in the univariate or multivariate time series
at vector of different lag values as well as they can be applied on an output object 
from a fitted model described in the description of the function \code{\link{BoxPierce}}.
}

\value{The Ljung and Box test statistic with the associated p-values 
for different lags based on the asymptotic chi-square distribution with \code{k^2(lags-fitdf)} degrees of freedom.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Ljung, G.M. and Box, G.E.P (1978). "On a Measure of Lack of Fit in Time Series Models". 
Biometrika, 65, 297-303.
}

\seealso{
\code{\link[stats]{Box.test}}, \code{\link{BoxPierce}}, \code{\link{MahdiMcLeod}}, 
\code{\link{Hosking}}, \code{\link{MahdiMcLeod}}, \code{\link{portest}}, \code{\link{GetResiduals}}.
}

\examples{
x <- rnorm(100)
LjungBox(x)                              ## univariate test
x <- cbind(rnorm(100),rnorm(100))
LjungBox(x)                              ## multivariate test      
##
##
## Annual flow of the river Nile at Aswan - 1871 to 1970
fit <- arima(Nile, c(1, 0, 1))
lags <- c(5, 10, 20)
## Apply the univariate test statistic on the fitted model 
LjungBox(fit, lags)            ## Correct (no need to specify fitdf) 
LjungBox(fit, lags, fitdf = 2) ## Correct 
## Apply the test statistic on the residuals and set fitdf = 2 
res <- resid(fit)
LjungBox(res, lags)             ## Wrong (fitdf is needed!)  
LjungBox(res, lags, fitdf = 2)  ## Correct 
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
LjungBox(fit,lags)                ## Correct (no need to specify fitdf)
## Apply the test statistic on the residuals where fitdf = 2
res <- ts((fit$resid)[-(1:2), ])
LjungBox(res,lags)                ## Wrong (fitdf is needed!)  
LjungBox(res,lags,fitdf = 2)      ## Correct 
##
##
## Monthly log stock returns of Intel corporation data: Test for ARCH Effects 
monthintel <- as.ts(monthintel)
LjungBox(monthintel)                         ## Usual test 
LjungBox(monthintel,sqrd.res=TRUE)  ## Test for ARCH effects
##
#### Write a function to fit a model: Apply portmanteau test on fitted obj with class "list"
## Example
FitModel <- function(data){
    fit <- ar.ols(data, intercept = TRUE, order.max = 2)
    fitdf <- 2
    res <- res <- ts((fit$resid)[-(1:2), ]) 
 list(res=res,fitdf=fitdf)
}
data(WestGerman)
DiffData <- matrix(numeric(3 * 91), ncol = 3)
  for (i in 1:3) 
    DiffData[, i] <- diff(log(WestGerman[, i]), lag = 1)
Fit <- FitModel(DiffData)
LjungBox(Fit) 
}

\keyword{Portmanteau Test}
\keyword{ ts}

