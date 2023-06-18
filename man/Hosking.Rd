\name{Hosking}
\alias{Hosking}

\title{The Modified Multivariate Portmanteau Test, Hosking (1980)}

\description{
The modified multivariate portmanteau test suggested by Hosking (1980).
}

\usage{
Hosking(obj,lags=seq(5,30,5),fitdf=0,sqrd.res=FALSE)
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
}

\value{
The multivariate test statistic suggested by Hosking (1980) and its associated p-values 
for different lags based on the asymptotic chi-square distribution with \code{k^2(lags-fitdf)} degrees of freedom.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Hosking, J. R. M. (1980). "The Multivariate Portmanteau Statistic". 
Journal of American Statistical Association, 75, 602-608.
}

\seealso{
\code{\link[stats]{Box.test}}, \code{\link{BoxPierce}}, \code{\link{LjungBox}}, \code{\link{MahdiMcLeod}}, 
\code{\link{LiMcLeod}}, \code{\link{portest}}, \code{\link{GetResiduals}}.
}

\examples{
x <- rnorm(100)
Hosking(x)                              ## univariate test
x <- cbind(rnorm(100),rnorm(100))
Hosking(x)                              ## multivariate test
##
##
## Quarterly, west German investment, income, and consumption from 1960 Q1  to 1982 Q4 
data(WestGerman)
DiffData <- matrix(numeric(3 * 91), ncol = 3)
  for (i in 1:3) 
    DiffData[, i] <- diff(log(WestGerman[, i]), lag = 1)
fit <- ar.ols(DiffData, intercept = TRUE, order.max = 2)
lags <- c(5,10)
## Apply the test statistic on the fitted model (fitdf will be automatically applied)
Hosking(fit,lags,fitdf = 2)                          ## Correct (no need to specify fitdf)
Hosking(fit,lags)                                    ## Correct
## Apply the test statistic on the residuals
res <- ts((fit$resid)[-(1:2), ])
Hosking(res,lags,fitdf = 2)                          ## Correct
Hosking(res,lags)                                    ## Wrong (fitdf is needed!)  
##
##
## Write a function to fit a model: Apply portmanteau test on fitted obj with class "list"
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
Hosking(Fit) 
}

\keyword{Portmanteau Test}


