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
              \code{"arima0"}, \code{"Arima"}, \code{"varest"}, \code{"FitAR"}, \code{"FitFGN"},
              \code{"garch"}, or \code{"fGARCH"}. 
              \code{obj} may also an object with class \code{"list"} (see details and following example).}
  \item{lags}{vector of lag auto-cross correlation coefficients used for \code{LiMcLeod} test.}
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
\code{\link{Hosking}}, \code{\link{gvtest}}, \code{\link{portest}}, \code{\link{GetResiduals}}, \code{\link[TSA]{tar}}
}

\examples{
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
LiMcLeod(fit,lags,order = 2)        ## Correct
LiMcLeod(fit,lags)                  ## Correct
## Apply the test statistic on the residuals
res <- ts((fit$resid)[-(1:2), ])
LiMcLeod(res,lags,order = 2)        ## Correct
LiMcLeod(res,lags)                  ## Wrong
##############################################################
## Write a function to fit a model 
## Apply portmanteau test on fitted obj with class "list"
##############################################################
## Example 1
FitModel <- function(data){
    fit <- ar.ols(data, intercept = TRUE, order.max = 2)
    order <- 2
    res <- res <- ts((fit$resid)[-(1:2), ]) 
 list(res=res,order=order)
}
Fit <- FitModel(DiffData)
LiMcLeod(Fit) 
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
LiMcLeod(Fit)
}
}
\keyword{htest}
\keyword{ ts}

