\name{MahdiMcLeod}
\alias{MahdiMcLeod}

\title{Generalized Variance Portmanteau Test}

\description{
New generalized variance portmanteau test based on the determinant of the Hosking's autocorrelation block 
Toeplitz matrix with order \eqn{m+1} given in the function \code{\link{ToeplitzBlock}}, 
where \eqn{m} represents the order of the block matrix.
Originally, the generalized variance portmanteau test, \code{MahdiMcLeod}, 
for univariate time series was derived by Pena and Rodriguez (2002)
based on the gamma distribution. 
Lin and McLeod (2006) proposed the Monte-Carlo version of this test and 
Mahdi and McLeod (2012) extended both methods to the multivariate case. 
Simulation results suggest that the Monte-Carlo version of 
\code{MahdiMcLeod} statistic is more accurate
and powerful than its competitors proposed by Box and Pierce (1970), Ljung and Box (1978),
and Pena and Rodriguez (2002, 2006) in the univariate time series and
Hosking (1980) and Li and McLeod (1981) in the multivariate time series.
}

\usage{
MahdiMcLeod(obj,lags=seq(5,30,5),order=0,season=1,squared.residuals=FALSE)
}

\arguments{
  \item{obj}{a univariate or multivariate series with class \code{"numeric"}, \code{"matrix"}, 
             \code{"ts"}, or \code{("mts" "ts")}.
             It can be also an object of fitted time-series model with class \code{"ar"}, 
              \code{"arima0"}, \code{"Arima"}, \code{("ARIMA" "Arima")}, \code{"lm"},
              \code{("glm" "lm")}, or \code{"varest"}. 
              \code{obj} may also an object with class \code{"list"} (see details and following examples).}
  \item{lags}{vector of lag auto-cross correlation coefficients used for \code{Hosking} test.}
  \item{order}{Default is zero for testing the randomness of a given sequence with 
              class \code{"numeric"}, \code{"matrix"}, \code{"ts"}, or \code{("mts" "ts")}. 
              In general \code{order} equals to the number of estimated parameters in the fitted model.
             If \code{obj} is an object with class \code{"ar"}, \code{"arima0"}, \code{"Arima"}, 
             \code{"varest"}, \code{("ARIMA" "Arima")}, or \code{"list"} then no need to enter 
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
\code{FitFGN()},  \code{fracdiff()}, \code{tar()}, etc.
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
Hosking, J. R. M. (1980). "The Multivariate Portmanteau Statistic". 
Journal of American Statistical Association, 75, 602-608.

Li, W. K. and McLeod, A. I. (1981). "Distribution of The Residual Autocorrelations in 
Multivariate ARMA Time Series Models". 
Journal of The Royal Statistical Society, Series B, 43, 231-239.

Lin, J.-W. and McLeod, A.I. (2006). "Improved Generalized Variance Portmanteau Test".
Computational Statistics and Data Analysis 51, 1731-1738.

Mahdi, E. and McLeod, A.I. (2012). "Improved Multivariate Portmanteau Test". 
Journal of Time Series Analysis, 33(2), 211-222.

Pena, D. and Rodriguez, J. (2002). "A Powerful Portmanteau Test of Lack of Test for Time Series".
Journal of American Statistical Association, 97, 601-610.

Pena, D. and Rodriguez, J. (2006). "The log of the determinant of the autocorrelation
matrix for testing goodness of fit in time series". Journal of Statistical Planning and
Inference, 136, 2706-2718.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link{ToeplitzBlock}}, \code{\link[stats]{Box.test}}, 
\code{\link{BoxPierce}}, \code{\link{LjungBox}}, \code{\link{Hosking}}, 
\code{\link{LiMcLeod}}, \code{\link{portest}}, \code{\link{GetResiduals}}.
}

\examples{
x <- rnorm(100)
MahdiMcLeod(x)                              ## univariate test
x <- cbind(rnorm(100),rnorm(100))
MahdiMcLeod(x)                              ## multivariate test      
##
##
## Annual flow of the river Nile at Aswan - 1871 to 1970
fit <- arima(Nile, c(1, 0, 1))
lags <- c(5, 10, 20)
## Apply the univariate test statistic on the fitted model 
MahdiMcLeod(fit, lags)            ## Correct (no need to specify order) 
MahdiMcLeod(fit, lags, order = 2) ## Correct 
## Apply the test statistic on the residuals and set order = 2 
res <- resid(fit)
MahdiMcLeod(res, lags)             ## Wrong (order is needed!)  
MahdiMcLeod(res, lags, order = 2)  ## Correct 
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
MahdiMcLeod(fit,lags)                ## Correct (no need to specify order)
## Apply the test statistic on the residuals where order = 2
res <- ts((fit$resid)[-(1:2), ])
MahdiMcLeod(res,lags)                ## Wrong (order is needed!)  
MahdiMcLeod(res,lags,order = 2)      ## Correct 
##
##
## Monthly log stock returns of Intel corporation data: Test for ARCH Effects 
monthintel <- as.ts(monthintel)
MahdiMcLeod(monthintel)                         ## Usual test 
MahdiMcLeod(monthintel,squared.residuals=TRUE)  ## Test for ARCH effects
##
##
## Test for seasonality
## Accidental Deaths in the US 1973 - 1978
seasonal.arima <- arima(USAccDeaths, order = c(0,1,1), seasonal = list(order = c(0,1,1)))
MahdiMcLeod(seasonal.arima, lags = 5, season = 12)
## Quarterly U.K. economic time series from 1957 Q3 to 1967 Q4
cd <- EconomicUK[,1]
cd.fit <- arima(cd,order=c(0,1,0),seasonal=list(order=c(0,1,1),period=4))
MahdiMcLeod(cd.fit, lags = c(5,10), season = 4)
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
MahdiMcLeod(Fit) 
detach(package:FitAR)
##
## Example 2
require("TSA")
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
MahdiMcLeod(Fit)
detach(package:TSA)
##
## Example 3
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
MahdiMcLeod(Fit) 
}

\keyword{Portmanteau Test}


