"portest" <-
function (obj,lags=seq(5,30,5),test=c("MahdiMcLeod","BoxPierce","LjungBox",
        "Hosking","LiMcLeod","other"),fn=NULL,squared.residuals=FALSE,MonteCarlo=TRUE,
         innov.dist=c("Gaussian","t","stable", "bootstrap"), ncores=1,nrep=1000,
        model=list(sim.model=NULL,fit.model=NULL),pkg.name=NULL,set.seed=123,...)
{
    if ((test == "other" && is.null(fn))|| (test != "other" && !is.null(fn)))
       stop("inconsistent match between fn and test")
    if (!is.null(fn)) 
       fn <- match.fun(fn) 
    dots <- list(...) 
    order <- dots$order               ## degrees of freedom of fitted model: e.g., p+q in ARMA models
    dft <- dots$dft                   ## df for t-distribution
    trunc.lag <- dots$trunc.lag       ## optional trunc lag pass varima.sim() function to trunc simulated series
    season <- dots$season             ## seasonal periodicity for testing seasonality
    if (is.null(order)) order <- 0
    if (is.null(season)) season <- 1
    innov.dist <- match.arg(innov.dist)
    test <- match.arg(test)
     TestType <- "0"
    if (class(obj) == "ts" || class(obj) == "numeric" || class(obj) ==
        "matrix" || (class(obj)[1] == "mts" && class(obj)[2] == "ts"))
        TestType <- "1"
    if (class(obj) == "ar" || class(obj) == "arima0" || class(obj) == 
        "Arima" || (class(obj)[1] == "ARIMA" && class(obj)[2] == "Arima") || class(obj) == "varest" || class(obj) == "lm"
        || (class(obj)[1] == "glm" && class(obj)[2] == "lm") || class(obj) == "list") 
        TestType<-"2"
    if (TestType == "0") 
        stop("obj must be class ar, arima0, Arima, (ARIMA Arima), varest, lm, (glm lm), ts, numeric, matrix, (mts ts), or list")
    if (TestType == "1") { ## apply the test on random series
        res <- as.ts(obj)
        Order <- order
    }
    else { ## apply the test on fitted models
        GetResid <- GetResiduals(obj)
        res <- GetResid$res
        Order <- GetResid$order
    }
    k <- NCOL(res)
    n <- NROW(res)
    res <- matrix(res,ncol=k,nrow=n)
    if (MonteCarlo == FALSE){ ## portmanteau test based on asymptotic distribution
      if (test =="MahdiMcLeod")
        return(MahdiMcLeod(res, lags, Order, season, squared.residuals))
      else if (test =="BoxPierce")
        return(BoxPierce(res, lags, Order, season, squared.residuals))
      else if (test =="LjungBox")
        return(LjungBox(res, lags, Order, season, squared.residuals))
      else if (test =="Hosking")
        return(Hosking(res, lags, Order, season, squared.residuals))
      else if (test =="LiMcLeod")
        return(LiMcLeod(res, lags, Order, season, squared.residuals))
      else if (test =="other")
        return(cbind(lags,fn(res, lags)))
   }
    else { ## portmanteau test based on Monte-Carlo procedures (compare observed.test with simulated.test)
          if (test =="MahdiMcLeod"){
              ans <- MahdiMcLeod(res, lags, Order, season, squared.residuals)
              obs.stat <- ans[, 2]
         }
         else if (test =="BoxPierce"){
             ans <- BoxPierce(res, lags, Order, season, squared.residuals)
             obs.stat <- ans[, 2]
        }
        else if (test =="LjungBox"){
            ans <- LjungBox(res, lags, Order, season, squared.residuals)
            obs.stat <- ans[, 2]
       }
       else if (test =="Hosking"){
           ans <- Hosking(res, lags, Order, season, squared.residuals)
           obs.stat <- ans[, 2]
      }
      else if (test =="LiMcLeod"){
           ans <- LiMcLeod(res, lags, Order, season, squared.residuals)
           obs.stat <- ans[, 2]
      }
      else if (test =="other"){
           ans <- cbind(lags,matrix(fn(res, lags),nrow=length(lags)),NA)
           obs.stat <- ans[,2]
     }
      if (TestType == 1)
           sigma <- matrix(stats::acf(res, lag.max = 1, plot = FALSE,type = "covariance")$acf[1, , ], k, k)
      else { ## we need to extract parameters from fitted model so we use them for simulation
                        if (all(class(obj) == "ar")) {##ar models on univariate series (Model 1) and multivariate series (Model2) 
                            p <- obj$order
                            q <- 0
                            if (k == 1)
                                Model <- 1
                            else
                                 Model <- 2
                            sigma <- obj$var.pred
                            if (is.array(obj$ar)) {
                                arrayphi <- array(numeric(k * k * p), dim = c(k^2, p))
                                for (i in 1:p) arrayphi[, i] <- c(obj$ar[i, ,])
                                phi <- array(c(arrayphi), dim = c(k, k, p))
                            }
                            else phi <- obj$ar
                            theta <- NULL
                            if (!is.null(obj$x.intercept)){
                                constant <- obj$x.intercept
                                intercept <- TRUE
                            }
                            else {
                                 constant <- rep(0, k)
                                 intercept <- FALSE
                            }
                            trend <- rep(0, k)
                            demean <- obj$x.mean
                            if(all(demean==0))
                               Demean <- FALSE
                            else
                               Demean <- TRUE
                            }
                else if (all(class(obj) == "varest")) {## multivariate var models - vars package
                              sigma <- summary(obj)[[3]]
                              p <- obj$p
                              q <- 0
                              theta <- NULL
                              if (obj$type == "none") {
                                  Model <- "3A"
                                  Phi <- matrix(numeric(p * k^2), nrow = k, ncol = p * k)
                                  constant <- NA
                                  trend <- NA
                                  demean <- NA
                                  for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                                  phi <- array(Phi, dim = c(k, k, p))
                              }
                              else if (obj$type == "const") {
                                  Model <- "3B"
                                  Phi <- matrix(numeric(k + p * k^2), nrow = k,ncol = p * k + 1)
                                  for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                                  constant <- Phi[, p * k + 1]
                                  trend <- NA
                                  demean <- NA
                                  phi <- array(Phi[, -(p * k + 1)], dim = c(k, k, p))
                              }
                              else if (obj$type == "trend") {
                                  Model <- "3C"
                                  Phi <- matrix(numeric(k + p * k^2), nrow = k,ncol = p * k + 1)
                                  for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                                  constant <- NA
                                  trend <- Phi[, p * k + 1]
                                  demean <- NA
                                  phi <- array(Phi[, -(p * k + 1)], dim = c(k, k, p))
                               }
                               else if (obj$type == "both") {
                                  Model <- "3D"
                                  Phi <- matrix(numeric(2 * k + p * k^2), nrow = k, ncol = p * k + 2)
                                  for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                                  constant <- Phi[, p * k + 1]
                                  trend <- Phi[, p * k + 2]
                                  demean <- NA
                                  phi <- array(Phi[, -((p * k + 1):(p * k + 2))], dim = c(k, k, p))
                                }
                         }
                else if (all(class(obj) == "arima0") || all(class(obj) == "Arima")) {## arima0 or arima models
                                Model <- 4
                                pdq <- obj$arma
                                if(as.integer(pdq[3])!= 0 || as.integer(pdq[4]) != 0 || as.integer(pdq[7]) != 0)
                                   stop ( "Use auto.arima() or Arima() function for seasonal arima models instead of arima() and arima0()")
	                            p <- pdq[1]
	                            q <- pdq[2]
	                            d <- pdq[6]
                                phi <- theta <- NULL  ## White Noise
                                if (p > 0 && q == 0)
                                    phi <- as.vector(obj$coef[1:p])
                                else if (p > 0 && q >0){
                                     phi <- as.vector(obj$coef[1:p])
                                     theta <- as.vector(obj$coef[(p+1):(p+q)])
                                }
                                else if (p==0 && q>0)
                                     theta <- as.vector(obj$coef[1:q])
                                sigma <- obj$sigma2
                                l <- length(obj$coef)
                                if (d==0){ ## arima with intercept
                                     demean <- as.vector(obj$coef[p+q+1])
                                     if (l > (p+q+1))
                                      stop ( "Use auto.arima() or Arima() function for arima models with xreg instead of arima() and arima0()")
                                }
                                else {
                                     demean <- 0
                                     if (l > (p+q))
                                       stop ( "Use auto.arima or Arima() function for arima models with xreg instead of arima() and arima0()")
                                }
                         }
                    else if (all (class(obj)[1] == "ARIMA" && class(obj)[2] == "Arima")){
                      Model <- 5
                      sigma <- obj$sigma2
                       pdq <- obj$arma
	                            p <- pdq[1]
	                            q <- pdq[2]
	                            ps <- pdq[3]
	                            qs <- pdq[4]
	                            season <- pdq[5]
	                            d <- pdq[6]
	                            ds <- pdq[7]
                               if (d==0){ 
                                   demean <- TRUE
                                   drift <- FALSE
                                   constant <- TRUE
                                }
                                else {
                                   demean <- FALSE
                                   drift <- TRUE
                                   constant <- FALSE
                                }               
                    }
                    else if (all(class(obj) == "lm") || all (class(obj)[1] == "glm" && class(obj)[2] == "lm")){
                        Model <- 6
                   }
                    else if (all(class(obj) == "list")) {
                           Model <- 7
                           stopifnot(!is.null(pkg.name)==TRUE)
                           pkg.name <<- as.name(pkg.name)
                           stopifnot(length(model)==2)
                       sim.model <- as.function(model[[1]])
                       fit.model <- as.function(model[[2]])
                     }
      }
      if (TestType == "1") { ## simulated test statistic applied on the residual 
          (OneMonteCarlo <- function() {
                if (innov.dist == "Gaussian") ## Gaussian distribution
                  rboot <- ts(matrix(stats::rnorm(k*n),ncol=k))
                else if (innov.dist == "t"){ ## t-distribution
                   if (is.null(dft)|| dft <= 0) 
                     stop("degrees of freedom for t-distribution is missing or entered as negative value") 
                   else
                     rboot <- ts(matrix(stats::rnorm(k*n),ncol=k)/sqrt(rchisq(k*n,dft)/dft))
                }
                else if (innov.dist == "stable"){ ## stable distribution
                   par.stable <- matrix(fitstable(res), ncol = 4)
            		Alpha<-par.stable[,1]
            		Beta<-par.stable[,2]
           		Scale<-par.stable[,3]
            		Location<-par.stable[,4]
                   rboot <- ts(rStable(n=n,Alpha=Alpha,Beta=Beta,Scale=Scale,Location=Location))
                }
                else if (innov.dist == "bootstrap") ## resampled residuals instead of distribution
                  rboot <- ts(matrix(sample(x=res, size=k*n, replace=TRUE),ncol = k))
              if (test =="MahdiMcLeod")
                  OneSim.stat <- MahdiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
              else if (test =="BoxPierce")
                  OneSim.stat <- BoxPierce(rboot, lags, Order, season, squared.residuals)[,2]
              else if (test =="LjungBox")
                  OneSim.stat <- LjungBox(rboot, lags, Order, season, squared.residuals)[,2]
              else if (test =="Hosking")
                  OneSim.stat <- Hosking(rboot, lags, Order, season, squared.residuals)[,2]
              else if (test =="LiMcLeod")
                  OneSim.stat <- LiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
              else if (test =="other")
                  OneSim.stat <-  matrix(fn(rboot, lags),nrow=length(lags))[,1]
              return(OneSim.stat)
          })
      }
      else { ## Simulated test statistic using the previously extracted fitted parameters Model 1 to Model 7
          (OneMonteCarlo <- function() { 
                if (Model == 1) {## univariate models with class ar
                if (innov.dist == "Gaussian")
                  innov <- function(n) ts(stats::rnorm(n, mean = demean, sd = sqrt(sigma)))
                else if (innov.dist == "t"){
                   if (is.null(dft)|| dft <= 0) 
                     stop("degrees of freedom for t-distribution is missing or entered as negative value") 
                  else 
                     innov <- function(n) ts(stats::rnorm(n, mean = demean, sd = sqrt(sigma))/sqrt(rchisq(n,dft)/dft))
                }
                else if (innov.dist == "stable"){ ## stable distribution
                   par.stable <- fitstable(res)
            		Alpha<-par.stable[1]
            		Beta<-par.stable[2]
           		      Scale<-par.stable[3]
            		Location<-par.stable[4]
                   innov <- function(n) ts(rStable(n=n,Alpha=Alpha,Beta=Beta,Scale=Scale,Location=Location))
                }
                else if (innov.dist == "bootstrap")
                  innov <- function(n) sample(x=res, size=n, replace=TRUE)
                 Sim.Data <- stats::arima.sim(n = n,list(ar = phi, ma = theta), rand.gen=innov)
                  FitSimModel <- stats::ar(Sim.Data, aic = FALSE, demean = Demean, order.max = p,method="yule-walker")
                  rboot <- FitSimModel$resid[-(1:p)]
                  if (test =="MahdiMcLeod")
                   OneSim.stat <- MahdiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="other")
                    OneSim.stat <-  matrix(fn(rboot, lags),nrow=length(lags))[,1]
                  return(OneSim.stat)
                }
                else if (Model == 2) {## multivariate models with class ar
			if (innov.dist == "Gaussian"){
				innov <- NULL
                        dft <- NULL
				par.stable <- NULL
				trunc.lag <- trunc.lag
			}
                else if (innov.dist == "t"){
                   if (is.null(dft)|| dft <= 0) 
                     stop("degrees of freedom for t-distribution is missing or entered as negative value") 
                   else 
                     innov <- NULL
                     dft <- dft
                     par.stable <- NULL
                     trunc.lag <- trunc.lag
                }
                else if (innov.dist == "stable"){ 
                      innov <- NULL
                      dft <- NULL
                      par.stable <- matrix(fitstable(res), ncol = 4)
                      trunc.lag <- trunc.lag
                }
                else if (innov.dist == "bootstrap"){
                      innov <- res
                      dft <- NULL
                      par.stable <- NULL
                      trunc.lag <- trunc.lag
                }
                    Sim.Data <- varima.sim(model=list(ar = phi, ma = theta), n = n, k=k, 
                       constant = constant, trend = trend, demean = demean, sigma = sigma,
                       innov = innov, innov.dist=innov.dist,dft=dft,par.stable=par.stable,trunc.lag=trunc.lag)
                  FitSimModel <- stats::ar.ols(Sim.Data, aic = FALSE, demean =Demean, intercept=intercept,order.max = p)
                  rboot <- ts(as.matrix(FitSimModel$resid)[-(1:p),])
                  if (test =="MahdiMcLeod")
                   OneSim.stat <- MahdiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="other")
                    OneSim.stat <- matrix(fn(rboot, lags),nrow=length(lags))[,1]
                  return(OneSim.stat)
                }
                else if (Model == "3A" || Model == "3B" || Model ==
                  "3C" || Model == "3D") {## multivariate models with class varest (vars package)
			if (innov.dist == "Gaussian"){
				innov <- NULL
                        dft <- NULL
				par.stable <- NULL
				trunc.lag <- trunc.lag
			}
                else if (innov.dist == "t"){
                   if (is.null(dft)|| dft <= 0) 
                     stop("degrees of freedom for t-distribution is missing or entered as negative value") 
                   else 
                     innov <- NULL
                     dft <- dft
                     par.stable <- NULL
                     trunc.lag <- trunc.lag
                }
                else if (innov.dist == "stable"){ ## stable distribution
                      innov <- NULL
                      dft <- NULL
                      par.stable <- matrix(fitstable(res), ncol = 4)
                      trunc.lag <- trunc.lag
                }
                else if (innov.dist == "bootstrap"){
                      innov <- res
                      dft <- NULL
                      par.stable <- NULL
                      trunc.lag <- trunc.lag
                }
                    Sim.Data <- varima.sim(model=list(ar = phi, ma = theta), n = n, k = k, 
                       constant = constant, trend = trend, demean = demean, sigma = sigma,
                       innov = innov, innov.dist=innov.dist,dft=dft,par.stable=par.stable,trunc.lag=trunc.lag)
                  if (Model == "3A")
                    FitSimModel <- vars::VAR(Sim.Data, p = p, type = "none")
                  else if (Model == "3B")
                    FitSimModel <- vars::VAR(Sim.Data, p = p, type = "const")
                  else if (Model == "3C")
                    FitSimModel <- vars::VAR(Sim.Data, p = p, type = "trend")
                  else if (Model == "3D")
                    FitSimModel <- vars::VAR(Sim.Data, p = p, type = "both")
                  rboot <- resid(FitSimModel)
                  if (test =="MahdiMcLeod")
                    OneSim.stat <- MahdiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="other")
                    OneSim.stat <- matrix(fn(rboot, lags),nrow=length(lags))[,1]
                  return(OneSim.stat)
                }
                else if (Model == 4) {## univariate models with class arima0 or Arima
                if (innov.dist == "Gaussian")
                  innov <- function(n) ts(stats::rnorm(n, mean = demean, sd = sqrt(sigma)))
                else if (innov.dist == "t"){
                   if (is.null(dft)|| dft <= 0) 
                     stop("degrees of freedom for t-distribution is missing or entered as negative value") 
                   else
                     innov <- function(n) ts(stats::rnorm(n, mean = demean, sd = sqrt(sigma))/sqrt(rchisq(n,dft)/dft))
                }
                else if (innov.dist == "stable"){ ## stable distribution
                   par.stable <- fitstable(res)
            		Alpha<-par.stable[1]
            		Beta<-par.stable[2]
           			Scale<-par.stable[3]
            		Location<-par.stable[4]
                   innov <- function(n) ts(rStable(n=n,Alpha=Alpha,Beta=Beta,Scale=Scale,Location=Location))
                }
                else if (innov.dist == "bootstrap")
                  innov <- function(n) sample(x=res, size=n, replace=TRUE)
                Sim.Data <- stats::arima.sim(n = n, list(ar = phi, ma = theta), rand.gen=innov)
                  if (d>0)
                   FitSimModel <- stats::arima(Sim.Data, order = c(p,d,q),include.mean=TRUE)
                  else
                   FitSimModel <- stats::arima(Sim.Data, order = c(p,d,q),include.mean=FALSE)
                  rboot <- FitSimModel$resid
                  if (test =="MahdiMcLeod")
                   OneSim.stat <- MahdiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="other")
                    OneSim.stat <- matrix(fn(rboot, lags),nrow=length(lags))[,1]
                 return(OneSim.stat)
                }
          else if (Model == 5) {## ## univariate models with class Arima arima -use simulate.Arima() function in forecast package
                if (innov.dist == "Gaussian")
                  innov <- ts(stats::rnorm(n, mean = demean, sd = sqrt(sigma)))
                else if (innov.dist == "t"){
                   if (is.null(dft)|| dft <= 0) 
                     stop("degrees of freedom for t-distribution is missing or entered as negative value") 
                   else
                     innov <- ts(stats::rnorm(n, mean = demean, sd = sqrt(sigma))/sqrt(rchisq(n,dft)/dft))
                }
                else if (innov.dist == "stable"){ 
                   par.stable <- fitstable(res)
            		Alpha<-par.stable[1]
            		Beta<-par.stable[2]
           		Scale<-par.stable[3]
            		Location<-par.stable[4]
                   innov <- ts(rStable(n=n,Alpha=Alpha,Beta=Beta,Scale=Scale,Location=Location))
                }
                else if (innov.dist == "bootstrap") 
                   innov <- sample(x=res, size=n, replace=TRUE)
                Sim.Data <- simulate(obj, nsim = length(obj$x),
                       seed = set.seed, xreg = obj$xreg, innov = innov, lambda = obj$lambda)
                       FitSimModel <- forecast::Arima(Sim.Data,order=c(p,d,q),seasonal=list(order=c(ps,ds,qs),period=season),xreg=obj$xreg,include.drift=drift,include.mean=demean,include.constant=constant)
                rboot <- FitSimModel$resid
                  if (test =="MahdiMcLeod")
                   OneSim.stat <- MahdiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="other")
                    OneSim.stat <- matrix(fn(rboot, lags),nrow=length(lags))[,1]
                  return(OneSim.stat)
          }
          else if (Model == 6) {## univariate models - lm() or glm() functions
              if (innov.dist == "Gaussian")
                  innov <- ts(stats::rnorm(n))
                else if (innov.dist == "t"){
                   if (is.null(dft)|| dft <= 0) 
                     stop("degrees of freedom for t-distribution is missing or entered as negative value") 
                   else
                     innov <- ts(stats::rnorm(n)/sqrt(rchisq(n,dft)/dft))
                }
                else if (innov.dist == "stable"){ 
                   par.stable <- fitstable(res)
            		Alpha<-par.stable[1]
            		Beta<-par.stable[2]
           			Scale<-par.stable[3]
            		Location<-par.stable[4]
                   innov <- ts(rStable(n=n,Alpha=Alpha,Beta=Beta,Scale=Scale,Location=Location))
                }
                else if (innov.dist == "bootstrap")
                  innov <- sample(x=res, size=n, replace=TRUE)
                 Sim.response <- ts(stats::simulate(obj, nsim=1))+innov
                 Sim.Data <- data.frame(obj$model[,-1],Sim.response)
                 colnames(Sim.Data) <- c(names(obj$model)[-1],"y")
               if (all(class(obj) == "lm"))
                    FitSimModel <- stats::lm(update.formula(formula(obj),y ~.), data=Sim.Data)
               else if (all(class(obj)[1] == "glm" && class(obj)[2] == "lm")) 
                FitSimModel <- stats::glm(update.formula(formula(obj),y ~.), data=Sim.Data)        
                 rboot <- FitSimModel$resid
                  if (test =="MahdiMcLeod")
                   OneSim.stat <- MahdiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="other")
                    OneSim.stat <- matrix(fn(rboot, lags),nrow=length(lags))[,1] 
                  return(OneSim.stat)
          }
          else if (Model == 7) {## univariate or multivariate models with class list
                  Sim.Data <- sim.model(obj)
                  FitSimModel <- fit.model(Sim.Data)
                  rboot <- ts(FitSimModel$res)
                  if (test =="MahdiMcLeod")
                   OneSim.stat <- MahdiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, season, squared.residuals)[,2]
                  else if (test =="other")
                    OneSim.stat <- matrix(fn(rboot, lags),nrow=length(lags))[,1]
                  return(OneSim.stat)
               }
            })
       }
     if(all(ncores == floor(ncores))==FALSE)
        stop ("ncores must be a positive integer >=1")
     if (ncores == 1) { ## Monte-Carlo with single CPU
          set.seed(set.seed)
         sim.stat <- replicate(nrep, OneMonteCarlo())
       }
     else{ ## Monte-Carlo with multiple CPUs 
           OneMonteCarlo <<- OneMonteCarlo
              if (all(class(obj) == "list")){
                sim.model <<- sim.model
                fit.model <<- fit.model
              }
            cl <- parallel::makeCluster(ncores)
              parallel::clusterSetRNGStream(cl, set.seed)
             if (all(class(obj) == "list")){
               parallel::clusterExport(cl, list("sim.model","fit.model"))
               Package <- as.call(list(library, as.name(pkg.name)))
               parallel::clusterCall(cl,eval,Package,env = .GlobalEnv)
             }
             parallel::clusterExport(cl, list("GetResiduals","MahdiMcLeod","LjungBox",
                "BoxPierce","Hosking","LiMcLeod","ImpulseVMA","InvertQ",
                "OneMonteCarlo", "varima.sim","vma.sim", "ToeplitzBlock", "simulate"))
             if (innov.dist == "stable"){
                parallel::clusterEvalQ(cl, library("akima"))
                parallel::clusterExport(cl, list("fitstable", "rStable"))
             }
            sim.stat <- parallel::parSapply(cl, 1:nrep, function(j) (OneMonteCarlo()))
            parallel::stopCluster(cl)
        }
       pvalue <- numeric(length(lags)) ## output for portmanteau tests
        for (i in 1:length(lags)) {
            if (is.matrix(sim.stat))
                pvalue[i] <- (1 + sum(as.numeric(sim.stat[i,] >= obs.stat[i])))/(nrep + 1)
            else 
                pvalue[i] <- (1 + sum(as.numeric(sim.stat >=obs.stat[i])))/(nrep + 1)
       }
       if (all(ans[,2]<0)) 
          ans2 <- -ans[,2]
       else 
          ans2 <- ans[,2]
       summary <- matrix(c(lags,ans2,ans[,3],pvalue),ncol=4)
       dimnames(summary) <- list(rep("", length(lags)),c("lags","statistic","df","p-value"))
      if (MonteCarlo == FALSE)
        return(summary)
      else
        return(summary[,-3])
     }
}
