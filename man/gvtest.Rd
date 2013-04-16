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
              \code{"arima0"}, \code{"Arima"}, \code{"varest"}, \code{"FitAR"}, \code{"FitFGN"},
              \code{"garch"}, or \code{"fGARCH"}. 
              \code{obj} may also an object with class \code{"list"} (see details and following example).}
  \item{lags}{vector of lag auto-cross correlation coefficients used for \code{gvtest} test.}
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
The generalized variance portmanteau test statistic and its associated p-values for different lags based on asymptotic 
chi-square as given in Mahdi and McLeod (2012).
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Mahdi, E. and McLeod, A.I. (2012). "Improved Multivariate Portmanteau Test". 
Journal of Time Series Analysis, 33(2), 211-222.

Pena, D. and Rodriguez, J. (2002). "A Powerful Portmanteau Test of Lack of Test for Time Series".
Journal of American Statistical Association, 97, 601-610.

Pena, D. and Rodriguez, J. (2006). "The log of the determinant of the autocorrelation
matrix for testing goodness of fit in time series". Journal of Statistical Planning and
Inference, 136, 2706-2718.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link[stats]{Box.test}}, \code{\link{BoxPierce}}, 
\code{\link{LjungBox}}, \code{\link{Hosking}}, 
\code{\link{LiMcLeod}}, \code{\link{portest}}, \code{\link{ToeplitzBlock}}, 
\code{\link{GetResiduals}}, \code{\link[TSA]{tar}}
}

\examples{
\dontrun{
x <- rnorm(100)
gvtest(x)            
##############################################################
## Measurements of the annual flow of the river Nile at Aswan 
## from the years 1871 to 1970:
##############################################################
fit <- arima(Nile, c(1, 0, 1))
lags <- c(5, 10, 20, 30)
## Apply the univariate test statistic on the fitted model 
gvtest(fit, lags)               ## Correct 
gvtest(fit, lags, order = 2)    ## Correct 
## Apply the test statistic on the residuals and set order = 2 
res <- resid(fit)
gvtest(res, lags)               ## Wrong  
gvtest(res, lags, order = 2)    ## Correct 
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
gvtest(fit,lags)                ## Correct 
## Apply the test statistic on the residuals where order = 2
res <- ts((fit$resid)[-(1:2), ])
gvtest(res,lags)                ## Wrong 
gvtest(res,lags,order = 2)      ## Correct 
##############################################################
## Monthly log stock returns of Intel corporation data
## Test for ARCH Effects 
##############################################################
monthintel <- as.ts(monthintel)
gvtest(monthintel)                ## Usual test
gvtest(monthintel,SquaredQ=TRUE)  ## Test for ARCH effects
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
gvtest(Fit) 
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
gvtest(Fit)
}
}

\keyword{htest}


