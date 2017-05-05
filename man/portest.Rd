\name{portest}
\alias{portest}

\title{Portmanteau Test Statistics}

\description{
 Univariate or multivariate portmanteau 
test statistics of \code{BoxPierce}, \code{MahdiMcLeod}, \code{Hosking}, 
\code{LiMcLeod}, \code{LjungBox}, and possibly any other test statistic using 
Monte-Carlo techniques or asymptotic distributions.
}

\usage{
 portest(obj,lags=seq(5,30,5),test=c("MahdiMcLeod","BoxPierce","LjungBox",
      "Hosking","LiMcLeod","other"),fn=NULL,squared.residuals=FALSE,MonteCarlo=TRUE,
         innov.dist=c("Gaussian","t","stable", "bootstrap"), ncores=1,nrep=1000,
        model=list(sim.model=NULL,fit.model=NULL),pkg.name=NULL,set.seed=123,...)
}

\arguments{
  \item{obj}{if \code{obj} is an object of class \code{"ar"}, \code{"arima0"}, \code{"Arima"}, 
             \code{("ARIMA" "Arima")}, \code{"lm"}, \code{("glm" "lm")}, 
             \code{"varest"}, or \code{"list"} 
             then a portmanteau goodness-of-fit test is done on the fitted model. 
             Otherwise, for \code{obj} with class \code{"ts"}, \code{"numeric"}, \code{"matrix"}, 
             or \code{("mts" "ts")}, a test of randomness is done.}
  \item{lags}{vector of lag values is used for portmanteau test.}
  \item{test}{portmanteau test statistic type.}
    \item{squared.residuals}{as described in \code{BoxPierce}, \code{MahdiMcLeod}, 
             \code{Hosking}, \code{LiMcLeod}, and \code{LjungBox}.} 
 \item{fn}{a function calculates the test statistic that is associated with \code{test = "other"}. 
   For example, \code{fn} can be a function returns the generalized Durbin-Watson test statistic values calculated at
different lags. This function has at least two inputs: \code{obj} and \code{lags}, where \code{obj}
and \code{lags} are described as above.} 
 \item{MonteCarlo}{if \code{TRUE} then apply the Monte-Carlo version of portmanteau statistic. 
             Otherwise, apply the asymptotic distribution. }
  \item{innov.dist}{distribution to generate univariate or multivariate innovation process. 
  This could be \code{Gaussian}, \code{t}, \code{stable}, 
  or \code{bootstrap} using resampled errors rather than distributed errors.
   Default is \code{Gaussian}.}
 \item{ncores}{number of cores needed to use in parallel calculations.
             Default is a single CPU.}
  \item{nrep}{number of replications needed for Monte-Carlo test.}
  \item{model}{additional argument defined as a list with two specified functions, 
             \code{sim.model} and \code{fit.model}.
             This argument is needed when the class of \code{obj} is \code{"list"} 
             (see details and following example).} 
  \item{pkg.name}{the name of the required library to be loaded if the Monte-Carlo significance test is used with
            an object \code{obj} with class \code{"list"}.} 
  \item{set.seed}{\code{set.seed} is initialized. Default seed is 123, but users can use any seed they wish.} 
  \item{...}{arguments to be passed to methods, such as \code{dft} degrees of freedom 
  needed to generate innovations with univariate/multivariate series with t-distribution innovations, 
 or \code{trunc.lag} used in \code{\link{varima.sim}} function, or \code{order} and \code{season} as described in
 \code{\link{BoxPierce}}, \code{\link{LjungBox}}, \code{\link{Hosking}}, \code{\link{LiMcLeod}} and \code{\link{MahdiMcLeod}}.}  
}

\details{
The portmanteau test statistics, \code{\link{MahdiMcLeod}}, \code{\link{BoxPierce}}, \code{\link{LjungBox}},
\code{\link{Hosking}}, and \code{\link{LiMcLeod}} are implemented based on the Monte-Carlo techniques
and the approximation asymptotic distributions as described in Mahdi and McLeod (2012). 
Any other possible test statistic is also implemented in this function by selecting the argument 
\code{test = "other"} and providing the test statistic as a function passing the argument \code{fn}.
The null hypothesis assuming that the fitted model is an adequate 
model and the residuals behave like white noise series. 
This function can be used for testing the  
adequacy in the nonseasonal fitted time series models.
this function can be used to check for randomness as well as to check for \code{ARCH}-\code{GARCH} effects.
Any other fitted model, for example, threashold autoregression model,
may also be tested for adequacy.
In this case, two functions, \code{sim.model()} and \code{fit.model()},
must be provided via the argument \code{func}.
The object \code{obj} is the output of the fitted model coded in
the function \code{fit.model} and it is a \code{"list"} with at least
\code{res}, the residuals from the fitted model in \code{fit.model()}, 
and \code{order}, the order of this fitted model.
The output from the function \code{sim.model()} is a simulated 
univariate or multivariate series from the fitted model obtained from the function
\code{fit.model()}.
The argument \code{pkg.name} represents the name of the \code{R} package where the fitted
model build in (see the last given example). 
The parallel computing using the \code{portes} package proposed by
Gonzalo Vera, Ritsert Jansen, and Remo Suppi (2008)
will run if one decide to choose the argument \code{MonteCarlo=TRUE} provided that 
\code{ncores} equals to a positive integer more than 1.
}

\value{
The portmanteau test statistic with the associated p-values 
for different lag values.
When the argument \code{MonteCarlo} is set to be \code{FALSE} then
the degrees of freedom will be an additional output.
}

\author{ Esam Mahdi and A.I. McLeod.}

\references{
Box, G.E.P. and Pierce, D.A. (1970). "Distribution of Residual Autocorrelation in Autoregressive-Integrated Moving Average
Time Series Models". Journal of American Statistical Association, 65, 1509-1526.

Chan KS, Ripley B (2012). TSA: Time Series Analysis. R package version 1.01, 
\url{https://CRAN.R-project.org/package=TSA}.

Fox, J and Weisberg, S and Adler, D and Bates, D and Baud-Bovy, G and Ellison, S and Firth, D and Friendly, M and 
Gorjanc, G and Graves, S and Heiberger, R and Laboissiere, R and Monette, G and Murdoch, D and Nilsson, H and Ogle, D 
and Ripley, B and Venables, W and Zeileis, A and R-Core (2016). car: Companion to Applied Regression. R package
version 2.1-4, \url{https://CRAN.R-project.org/package=car}.

Fraley C, Leisch F, Maechler M, Reisen V, Lemonte A (2012). fracdiff: Fractionally differenced ARIMA aka ARFIMA(p,d,q) models. 
R package version 1.4-2, \url{https://CRAN.R-project.org/package=fracdiff}.

Hosking, J. R. M. (1980). "The Multivariate Portmanteau Statistic". 
Journal of American Statistical Association, 75, 602-608.

John Haslett and Adrian E. Raftery (1989). "Space-time Modelling with Long-memory Dependence: 
Assessing Ireland's Wind Power Resource (with Discussion)". Applied Statistics, 38, 1-50.

Li, W. K. and McLeod, A. I. (1981). "Distribution of The Residual Autocorrelations in 
Multivariate ARMA Time Series Models". Journal of The Royal Statistical Society, Series B, 43, 231-239.

Lin, J.-W. and McLeod, A.I. (2006). "Improved Generalized Variance Portmanteau Test".
Computational Statistics and Data Analysis 51, 1731-1738.

Lin, J.-W. and McLeod, A.I. (2008). "Portmanteau Tests for ARMA Models with Infinite Variance".  
Journal of Time Series Analysis, 29, 600-617.

Ljung, G.M. and Box, G.E.P (1978). "On a Measure of Lack of Fit in Time Series Models". 
Biometrika, 65, 297-303.

Mahdi, E. and McLeod, A.I. (2012). "Improved Multivariate Portmanteau Test". 
Journal of Time Series Analysis, 33(2), 211-222.

McLeod A.I, Li W.K (1983). "Distribution of the Residual Autocorrelation in Multivariate
ARMA Time Series Models". Journal of Time Series Analysis, 4, 269-273.

McLeod A, Veenstra J (2014). FGN: Fractional Gaussian Noise and power law decay time series 
model fitting. R package version 2.0-12, \url{https://CRAN.R-project.org/package=FGN}.

McLeod A, Zhang Y, Xu C (2013). FitAR: Subset AR Model Fitting. R package version 1.94,
\url{https://CRAN.R-project.org/package=FitAR}.

Pena, D. and Rodriguez, J. (2002). "A Powerful Portmanteau Test of Lack of Test for Time Series".
Journal of American Statistical Association, 97, 601-610.

Pena, D. and Rodriguez, J. (2006). "The log of the determinant of the autocorrelation
matrix for testing goodness of fit in time series". Journal of Statistical Planning and
Inference, 136, 2706-2718.

Pfaff B, Stigler M (2013). vars: VAR Modelling. R package version 1.5-2, 
\url{https://CRAN.R-project.org/package=vars}.

Rob J Hyndman with contributions from George Athanasopoulos Slava Razbash DSZZYKCB,
Wang E (2017). forecast: Forecasting Functions for Time Series and Linear Models. 
R package version 8.0, \url{https://CRAN.R-project.org/package=forecast}.

Tierney, L., Rossini, A. J., Li, N., and Sevcikova, H. (2016). snow: Simple Network of Workstations.
\code{R} package version 0.4-2. \url{https://CRAN.R-project.org/package=snow}.

Trapletti A, Hornik K, LeBaron B (2017). tseries: Time Series Analysis and Computational Finance. 
R package version 0.10-38, \url{https://CRAN.R-project.org/package=tseries}.

Gonzalo Vera and Ritsert C. Jansen and Remo L. Suppi (2008). R/parallel - speeding 
up bioinformatics analysis with R. BMC Bioinformatics, 9:390.

Wuertz D, core team members R (2016). fGarch: Rmetrics - Autoregressive Conditional
Heteroskedastic Modelling. R package version 3010.82.1, \url{https://CRAN.R-project.org/package=fGarch}.
}

\seealso{
\code{\link[stats]{acf}}, \code{\link{ar}}, \code{\link{ar.ols}}, \code{\link{ar.burg}}, 
\code{\link{ar.yw}}, \code{\link{ar.mle}}, \code{\link{arima0}}, \code{\link{arima}}, 
\code{\link[stats]{lm}}, \code{\link[stats]{glm}}, 
\code{\link[stats]{Box.test}}, \code{\link{BoxPierce}}, \code{\link{LjungBox}}, \code{\link{MahdiMcLeod}}, 
\code{\link{LiMcLeod}}, \code{\link{portest}}, \code{\link{ToeplitzBlock}}, \code{\link{GetResiduals}}, 
\code{\link[TSA]{tar}}, \code{\link[forecast]{Arima}}, \code{\link[forecast]{auto.arima}}, 
\code{\link[vars]{VAR}}, \code{\link[FitAR]{FitAR}}, \code{\link[FGN]{FitFGN}}, \code{\link[fracdiff]{fracdiff}}, 
\code{\link[fGarch]{garchFit}}, \code{\link[tseries]{garch}}, \code{\link{varima.sim}}, \code{\link{fitstable}}.
}

\examples{
\dontrun{
#################################################################################
####                                                                         ####
####                             Portmanteau Tests                           ####  
####                                                                         #### 
#################################################################################
## Monte-Carlo (MC) and asymptotic tests for randomness series                 ##
#################################################################################
data("DEXCAUS")
returns <- log(DEXCAUS[-1]/DEXCAUS[-length(DEXCAUS)])
portest(returns)            ## MC using one CPU takes about 24.16 seconds
portest(returns, ncores=4)    ## MC using 4 CPUs takes about 9.51 seconds
portest(returns, MonteCarlo=FALSE)              ## asymptotic MahdiMcLeod
portest(returns,test="LjungBox", MonteCarlo=FALSE) ## asymptotic LjungBox 
#################################################################################
## Monte-Carlo goodness-of-fit arima test using 4 CPUs                         ##
#################################################################################
## arima() function takes about 11.32 seconds
## Example 1
ans1 <- arima(WWWusage,order=c(3,1,0))
portest(ans1, ncores = 4)
#
## arima0() function takes about 11.1 seconds
## Example 2
ans2 <- arima0(WWWusage,order=c(3,1,0))
portest(ans2, ncores = 4)
#
## Arima() or auto.Arima() function takes about 12.1 seconds
## Example 3
require("forecast")
ans3 <- Arima(WWWusage,order=c(3,1,0))
portest(ans3, ncores = 4)
detach(package:forecast)
#
## ar() function takes about 7.39 seconds
## Example 4
ans4 <- ar(Nile,order.max=2)
portest(ans4, ncores = 4)
#
## ar() function with your own R code takes about 8.75 seconds
## Example 5
fit.model <- function(data){
    fit <- ar(data,aic = FALSE, order.max=2)
    order <- 2
    res <- ts(fit$resid[-(1:order)]) 
    phi <- fit$ar
    theta <- NULL
    sigma <- fit$var.pred
    demean <- fit$x.mean
 list(res=res,phi=phi,theta=theta,sigma=sigma,demean=demean)
}
sim.model <- function(parSpec){
    res <- parSpec$res
    n <- length(res)  
    innov <- function(n) ts(stats::rnorm(n, mean = demean, sd = sqrt(sigma))) 
    phi <- parSpec$phi
    theta <- parSpec$theta
    sigma <- parSpec$sigma 
    demean <- parSpec$demean 
   arima.sim(n = n, list(ar = phi, ma = theta), rand.gen=innov)
}
ans5 <- fit.model(Nile)
portest(ans5,ncores=4,model=list(sim.model=sim.model,fit.model=fit.model),pkg.name="stats")
#################################################################################
## Monte-Carlo test for seasonality                                            ##
#################################################################################
## Accidental Deaths in the US 1973 - 1978
require("forecast")
seasonal.arima<-Arima(USAccDeaths,order=c(0,1,1),seasonal=list(order= c(0,1,1)))
portest(seasonal.arima,ncores=4,nrep=1000,lags=1:5)
## Quarterly U.K. economic time series from 1957 Q3 to 1967 Q4
cd <- EconomicUK[,1]
cd.fit <- Arima(cd,order=c(0,1,0),seasonal=list(order=c(0,1,1),period=4))
portest(cd.fit, lags = c(5,10),ncores=4)
detach(package:forecast)
#################################################################################
## Monte-Carlo test for linear models and time series regression               ##
#################################################################################
## Linear Model
require("car")
fit <- lm(fconvict ~ tfr + partic + degrees + mconvict, data=Hartnagel)
portest(fit,lags=1:3,ncores=4)          ## MC of MahdiMcLeod test
## MC of generalized Durbin-Watson test needs the argument function fn() as follows
fn <- function(obj,lags){
     test.stat <- numeric(length(lags))
       for (i in 1:length(lags))
          test.stat[i] <- -sum(diff(obj,lag=lags[i])^2)/sum(obj^2)
       test.stat
}	
portest(fit,lags=1:3,test="other",fn=fn,ncores=4) 
detach(package:car)
## Time series regression 
require("forecast")
fit.arima <- Arima(LakeHuron, order = c(2,0,0), xreg = time(LakeHuron)-1920)
portest(fit.arima,ncores=4)
detach(package:forecast)
#################################################################################
## Monte-Carlo goodness-of-fit VAR test - Multivariate series                  ##
#################################################################################
data("IbmSp500")
ibm <- log(IbmSp500[,2]+1)*100
sp500 <- log(IbmSp500[,3]+1)*100
IBMSP500 <- data.frame(cbind(ibm,sp500))
## ar.ols() function takes about 9.11 seconds
ans6 <-  ar.ols(IBMSP500, aic=FALSE, intercept=TRUE, order.max=5)
portest(ans6,nrep=100,test="MahdiMcLeod",ncores=4,innov.dist="t",dft=5)
## VAR() function takes about 11.55 seconds
require("vars")
ans7 <- VAR(IBMSP500, p=5) 
portest(ans7,nrep=100,test="MahdiMcLeod",ncores=4,innov.dist="bootstrap")
portest(ans7,test="Hosking",MonteCarlo=FALSE) ## asymptotic Hosking test
detach(package:vars)
#################################################################################
## Monte-Carlo test for models with heavy tails stable distributions           ##
#################################################################################
## It takes about 32.7 seconds on personal PC with 4 CPUs 
data("CRSP")
CRSP.AR5<- arima(CRSP, c(5, 0, 0))
lags <- c(10, 20, 30)
portest(CRSP.AR5,lags=lags,ncores=4,nrep=1000,innov.dist="stable")   
#################################################################################
## Monte-Carlo test for ARCH/GARCH effects using 4 CPUs                        ##
#################################################################################
## It takes about 12.65 seconds
data("monthintel")
returns <- as.ts(monthintel)
lags <- c(5, 10, 20, 40)
portest(returns, lags = lags, ncores = 4, squared.residuals = FALSE)
portest(returns,lags=lags,ncores=4,squared.residuals=TRUE,innov.dist="t",dft=5)
#################################################################################
## Monte-Carlo test for Threshold Autoregressive (TAR) Models                  ##
## It takes about 54.27 seconds on personal PC with 4 CPUs                     ##
#################################################################################
require("TSA")
fit.model <- function(data){
    fit <- TSA::tar(y=log(data),p1=4,p2=4,d=3,a=0.1,b=0.9,print=FALSE)
    res <- ts(fit$std.res)
    parSpec <- list(res=res,fit=fit)
  parSpec
}
sim.model <- function(parSpec){
    fit <- parSpec$fit   
  exp(tar.sim(fit)$y)
}
data(prey.eq)
portest(fit.model(prey.eq),ncores=4,model=list(sim.model,fit.model),pkg.name="TSA")
detach(package:TSA)
}
}

\keyword{Portmanteau Test}



