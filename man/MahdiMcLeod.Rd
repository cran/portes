\name{MahdiMcLeod}
\alias{MahdiMcLeod}

\title{Generalized Variance Portmanteau Test}

\description{
New generalized variance portmanteau test based on the determinant of the Hosking's autocorrelation block 
Toeplitz matrix with fitdf \eqn{m+1} given in the function \code{\link{ToeplitzBlock}}, 
where \eqn{m} represents the fitdf of the block matrix.
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
MahdiMcLeod(obj,lags=seq(5,30,5),fitdf=0,sqrd.res=FALSE)
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
MahdiMcLeod(fit, lags)            ## Correct (no need to specify fitdf) 
MahdiMcLeod(fit, lags, fitdf = 2) ## Correct 
## Apply the test statistic on the residuals and set fitdf = 2 
res <- resid(fit)
MahdiMcLeod(res, lags)             ## Wrong (fitdf is needed!)  
MahdiMcLeod(res, lags, fitdf = 2)  ## Correct 
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
MahdiMcLeod(fit,lags)                ## Correct (no need to specify fitdf)
## Apply the test statistic on the residuals where fitdf = 2
res <- ts((fit$resid)[-(1:2), ])
MahdiMcLeod(res,lags)                ## Wrong (fitdf is needed!)  
MahdiMcLeod(res,lags,fitdf = 2)      ## Correct 
##
##
## Monthly log stock returns of Intel corporation data: Test for ARCH Effects 
monthintel <- as.ts(monthintel)
MahdiMcLeod(monthintel)                         ## Usual test 
MahdiMcLeod(monthintel,sqrd.res=TRUE)  ## Test for ARCH effects
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
MahdiMcLeod(Fit) 
}

\keyword{Portmanteau Test}


