"portest" <-
function (obj,lags=seq(5,30,5),order=0,test=c("gvtest","BoxPierce","LjungBox",
        "Hosking","LiMcLeod"),MonteCarlo=TRUE,nslaves=1,NREP=1000,
        InfiniteVarianceQ=FALSE,SquaredQ=FALSE,SetSeed=TRUE) 
{
    test <- match.arg(test)
     TestType <- "0"
    if (class(obj) == "ts" || class(obj) == "numeric" || class(obj) == 
        "matrix" || (class(obj)[1] == "mts" && class(obj)[2] == 
        "ts")) 
        TestType <- "1"
    if (class(obj) == "ar" || class(obj) == "arima0" || class(obj) == 
        "Arima" || class(obj) == "varest" || class(obj) == "FitAR" || 
        class(obj) == "FitFGN") 
        TestType <- "2"
    if (TestType == "0") 
        stop("obj must be class ar, arima0, Arima, varest, FitAR, 
             FitFGN, ts, numeric, matrix, or (mts ts)")
    if (TestType == "1") {
        res <- as.ts(obj)
        Order <- order
    }
    else {
        GetResid <- GetResiduals(obj)
        res <- GetResid$res
        Order <- GetResid$order
    }
    k <- NCOL(res)
    n <- NROW(res)
    if (MonteCarlo == FALSE){ 
      if (test =="gvtest")
        return(gvtest(res, lags, Order, SquaredQ))
      else if (test =="BoxPierce")
        return(BoxPierce(res, lags, Order, SquaredQ))  
      else if (test =="LjungBox")
        return(LjungBox(res, lags, Order, SquaredQ)) 
      else if (test =="Hosking")
        return(Hosking(res, lags, Order, SquaredQ)) 
      else if (test =="LiMcLeod")
        return(LiMcLeod(res, lags, Order, SquaredQ)) 
   }           
    else {
        Trunc.Series <- min(100, ceiling(n/3))
        if (test =="gvtest"){
          ans <- gvtest(res, lags, Order, SquaredQ)
          obs.stat <- ans[, 2]
        }
        else if (test =="BoxPierce"){
          ans <- BoxPierce(res, lags, Order, SquaredQ)
          obs.stat <- ans[, 2]
        }
        else if (test =="LjungBox"){
          ans <- LjungBox(res, lags, Order, SquaredQ)
          obs.stat <- ans[, 2]
        }
        else if (test =="Hosking"){
          ans <- Hosking(res, lags, Order, SquaredQ)
          obs.stat <- ans[, 2]
        }
        else if (test =="LiMcLeod"){
          ans <- LiMcLeod(res, lags, Order, SquaredQ)
          obs.stat <- ans[, 2]
        }
      if (all(class(obj) == "ar")) {
            p <- obj$order
            q <- 0
            if (k == 1) 
                Model <- 1
            else Model <- 2
            sigma <- obj$var.pred
            if (is.array(obj$ar)) {
                arrayphi <- array(numeric(k * k * p), dim = c(k^2, 
                  p))
                for (i in 1:p) arrayphi[, i] <- c(obj$ar[i, , 
                  ])
                phi <- array(c(arrayphi), dim = c(k, k, p))
            }
            else phi <- obj$ar
            theta <- NULL
            if (!is.null(obj$x.intercept)) 
                constant <- obj$x.intercept
            else constant <- rep(0, k)
            trend <- rep(0, k)
            demean <- obj$x.mean
        }
        else if (all(class(obj) == "varest")) {
            sigma <- summary(obj)[[3]]
            p <- obj$p
            q <- 0
            theta <- NULL
            if (obj$type == "none") {
                Model <- "3A"
                Phi <- matrix(numeric(p * k^2), nrow = k, ncol = p * 
                  k)
                constant <- NA
                trend <- NA
                demean <- NA
                for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                phi <- array(Phi, dim = c(k, k, p))
            }
            else if (obj$type == "const") {
                Model <- "3B"
                Phi <- matrix(numeric(k + p * k^2), nrow = k, 
                  ncol = p * k + 1)
                for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                constant <- Phi[, p * k + 1]
                trend <- NA
                demean <- NA
                phi <- array(Phi[, -(p * k + 1)], dim = c(k, 
                  k, p))
            }
            else if (obj$type == "trend") {
                Model <- "3C"
                Phi <- matrix(numeric(k + p * k^2), nrow = k, 
                  ncol = p * k + 1)
                for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                constant <- NA
                trend <- Phi[, p * k + 1]
                demean <- NA
                phi <- array(Phi[, -(p * k + 1)], dim = c(k, 
                  k, p))
            }
            else if (obj$type == "both") {
                Model <- "3D"
                Phi <- matrix(numeric(2 * k + p * k^2), nrow = k, 
                  ncol = p * k + 2)
                for (i in 1:k) Phi[i, ] <- coef(obj)[[i]][, 1]
                constant <- Phi[, p * k + 1]
                trend <- Phi[, p * k + 2]
                demean <- NA
                phi <- array(Phi[, -((p * k + 1):(p * k + 2))], 
                  dim = c(k, k, p))
            }
        }
       else if (all(class(obj) == "arima0") || all(class(obj) == 
            "Arima")) {
          Model <- 4
          pdq <- obj$arma
          if(as.integer(pdq[3])!= 0 || as.integer(pdq[4]) != 0 || as.integer(pdq[7]) != 0) 
           stop("portes is applied only for nonseasonal arima models")
	    p <- pdq[1]
	    q <- pdq[2]
	    d <- pdq[6]
           phi <- theta <- NULL
            if (p > 0 && q == 0)               
                     phi <- as.vector(obj$coef[1:p])
            else if (p > 0 && q >0){
                     phi <- as.vector(obj$coef[1:p])
                     theta <- as.vector(obj$coef[(p+1):(p+q)])
            }
            else if (p==0 && q>0)
                     theta <- as.vector(obj$coef[1:q])
            sigma <- obj$sigma2
                   if (length(obj$coef)==p+q+2) 
                    constant <- as.vector(obj$coef[p+q+2])
                   else
                    constant <- 0  
                    trend <- 0
                   if (d==0) 
                     demean <- as.vector(obj$coef[p+q+1])
                   else
                     demean <- 0   
       }
       else if (all(class(obj) == "FitAR")) {
            Model <- 5
            p <- length(obj$phiHat)
            q <- 0
            phi <- obj$phiHat
            theta <- NULL
            sigma <- obj$sigsqHat
            constant <- 0
            trend <- 0
            demean <- obj$muHat
       }
       else if (all(class(obj) == "FitFGN")) {
            Model <- 6
            phi <- theta <- NULL
            p <- 1
            q <- 0
            H <- obj$H
            sigma <- sqrt(obj$sigsq)
            constant <- 0
            trend <- 0
            demean <- obj$muHat
       }
       else sigma <- matrix(acf(res, lag.max = 1, plot = FALSE, 
            type = "covariance")$acf[1, , ], k, k)
        if (InfiniteVarianceQ) {
            StableParameters <- matrix(fitstable(res), ncol = 4)
            ALPHA <- StableParameters[, 1]
            BETA <- StableParameters[, 2]
            GAMMA <- StableParameters[, 3]
            DELTA <- StableParameters[, 4]
            ALPHA <<- ALPHA
            BETA <<- BETA
            GAMMA <<- GAMMA
            DELTA <<- DELTA
       }
       else StableParameters <- NA
        if (TestType == "1") {
            (OneMonteCarlo <- function() {
                if (InfiniteVarianceQ) 
                  rboot <- ts(rstable(n, ALPHA, BETA, GAMMA, 
                    DELTA))
                else rboot <- t(crossprod(chol(sigma), matrix(rnorm(k * 
                  n), ncol = n)))
                if (test =="gvtest")
                  OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[, 
                  2]
                else if (test =="BoxPierce")
                  OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[, 
                  2]
                else if (test =="LjungBox")
                  OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[, 
                  2]
                else if (test =="Hosking")
                  OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[, 
                  2]
                else if (test =="LiMcLeod")
                  OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[, 
                  2]
                return(OneSim.stat)
            })
       }
       else {
            (OneMonteCarlo <- function() {
                if (Model == 1) {
                  if (InfiniteVarianceQ) 
                    Sim.Data <- varima.sim(phi = phi, theta = theta, d=0,
                      sigma = sigma, n = n, constant = constant, 
                      trend = trend, demean = demean, StableParameters = StableParameters, 
                      Trunc.Series = Trunc.Series)
                  else Sim.Data <- constant + arima.sim(n = n, 
                    list(ar = phi, ma = theta), sd = sqrt(sigma), 
                    mean = demean)
                  FitSimModel <- ar(Sim.Data, aic = FALSE, order.max = p)
                  rboot <- FitSimModel$resid[-(1:p)]
                  if (test =="gvtest")
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[, 
                   2]
                  return(OneSim.stat)
                }
                else if (Model == 2) {
                  Sim.Data <- varima.sim(phi = phi, theta = theta,d=NA, 
                    sigma = sigma, n = n, constant = constant, 
                    trend = trend, demean = demean, StableParameters = StableParameters, 
                    Trunc.Series = Trunc.Series)
                  FitSimModel <- ar.ols(Sim.Data, aic = FALSE, 
                    order.max = p)
                  rboot <- ts(as.matrix(FitSimModel$resid)[-(1:p), 
                    ])
                  if (test =="gvtest")
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[, 
                   2]
                  return(OneSim.stat)
                }
                else if (Model == "3A" || Model == "3B" || Model == 
                  "3C" || Model == "3D") {
                  Sim.Data <- varima.sim(phi = phi, theta = theta, d=NA,
                    sigma = sigma, n = n, constant = constant, 
                    trend = trend, demean = demean, StableParameters = StableParameters, 
                    Trunc.Series = Trunc.Series)
                  if (Model == "3A") 
                    FitSimModel <- VAR(Sim.Data, p = p, type = "none")
                  else if (Model == "3B") 
                    FitSimModel <- VAR(Sim.Data, p = p, type = "const")
                  else if (Model == "3C") 
                    FitSimModel <- VAR(Sim.Data, p = p, type = "trend")
                  else if (Model == "3D") 
                    FitSimModel <- VAR(Sim.Data, p = p, type = "both")
                  rboot <- resid(FitSimModel)
                  if (test =="gvtest")
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[, 
                   2]
                  return(OneSim.stat)
                }
                else if (Model == 4) {
                  if (InfiniteVarianceQ) 
                    Sim.Data <- varima.sim(phi = phi, theta = theta, d=d,
                      sigma = sigma, n = n, constant = constant, 
                      trend = trend, demean = demean, StableParameters = StableParameters, 
                      Trunc.Series = Trunc.Series)
                  else 
                      Sim.Data <- constant + arima.sim(n = n, list(ar = phi, ma = theta), sd = sqrt(sigma), mean = demean)
                  if (d>0)
                      Sim.Data <- diffinv(Sim.Data,differences = d)[-(1:d)]
                  if (constant !=0)
                      FitSimModel <- Arima(Sim.Data, order = c(p,d,q),include.drift = TRUE)                       
                  else                             
                     FitSimModel <- arima(Sim.Data, order = c(p,d,q))
                  rboot <- FitSimModel$resid
                  if (test =="gvtest")
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[, 
                   2]
                  return(OneSim.stat)
                }
                else if (Model == 5) {
                  if (InfiniteVarianceQ) 
                    Sim.Data <- varima.sim(phi = phi, theta = theta,d=0, 
                      sigma = sigma, n = n, constant = constant, 
                      trend = trend, demean = demean, StableParameters = StableParameters, 
                      Trunc.Series = Trunc.Series)
                  else Sim.Data <- arima.sim(n = n, list(ar = phi, 
                    ma = theta), sd = sqrt(sigma), mean = demean)
                  pvec <- obj$pvec
                  if (obj$SubsetQ) {
                    if (obj$ARModel=="ARz")
                        rboot<-ts(GetFitARz(Sim.Data, pvec=pvec, MeanValue=mean(Sim.Data))$res)
                    else
                        rboot<-ts(GetFitARpLS(Sim.Data, pvec=pvec)$res) #let function do mean-correction
                        }
                  else #must be AR(p)
                    rboot<-ts(GetFitARpLS(Sim.Data, pvec=pvec)$res)
                  if (test =="gvtest")
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[, 
                   2]
                  return(OneSim.stat)
                }
                else if (Model == 6) {
                  Sim.Data <- demean + SimulateFGN(n, H) * sigma
                  FitSimModel <- FitFGN(Sim.Data)
                  rboot <- ts(FitSimModel$res)
                  if (test =="gvtest")
                   OneSim.stat <- gvtest(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="BoxPierce")
                    OneSim.stat <- BoxPierce(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LjungBox")
                    OneSim.stat <- LjungBox(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="Hosking")
                    OneSim.stat <- Hosking(rboot, lags, Order, SquaredQ)[, 
                    2]
                  else if (test =="LiMcLeod")
                    OneSim.stat <- LiMcLeod(rboot, lags, Order, SquaredQ)[, 
                   2]
                  return(OneSim.stat)
                }
            })
       }
       if (nslaves == 1) {
        if (SetSeed) 
          set.seed(21597341)
         sim.stat <- replicate(NREP, OneMonteCarlo())
       }
       else {
            n <<- n
            k <<- k
            Order <<- Order
            obj <<- obj
            lags <<- lags
            res <<- res
            NREP <<- NREP
            if (TestType == "2") {
                Model <<- Model
                p <<- p
                phi <<- phi
                theta <<- theta
                constant <<- constant
                trend <<- trend
                demean <<- demean
            }
            sigma <<- sigma
            SquaredQ <<- SquaredQ
            InfiniteVarianceQ <<- InfiniteVarianceQ
            Trunc.Series <<- Trunc.Series
            OneMonteCarlo <<- OneMonteCarlo
            cl <- makeCluster(nslaves)
            if (SetSeed) 
              clusterSetupRNG(cl, seed = 21597341)
            clusterEvalQ(cl, library("akima"))
            clusterEvalQ(cl, library("FitAR"))
            if (all(class(obj) == "Arima")) 
                clusterEvalQ(cl, library("forecast"))
            if (all(class(obj) == "varest")) 
                clusterEvalQ(cl, library("vars"))
            if (all(class(obj) == "FitFGN")) 
                clusterEvalQ(cl, library("FGN"))
            clusterExport(cl, list("k", "n", "obj", "Order", 
                "sigma", "lags", "NREP", "res", "InfiniteVarianceQ", 
                "SquaredQ", "Trunc.Series"))
            clusterExport(cl, list("GetResiduals","gvtest","LjungBox",
                "BoxPierce","Hosking","LiMcLeod","ImpulseVMA","InvertQ", 
                "OneMonteCarlo", "varima.sim","vma.sim","ToeplitzBlock"))
            if (InfiniteVarianceQ) 
                clusterExport(cl, list("ALPHA", "BETA", "GAMMA", 
                  "DELTA", "interpp", "interpp.old", "fitstable", 
                  "rstable"))
            if (TestType == "2") 
                clusterExport(cl, list("p", "Model", "phi", "theta", 
                  "constant", "trend", "demean"))
            sim.stat <- parSapply(cl, 1:NREP, function(j) (OneMonteCarlo()))
            stopCluster(cl)
       }
       pvalue <- numeric(length(lags))
        for (i in 1:length(lags)) {
            if (is.matrix(sim.stat)) 
                pvalue[i] <- (1 + sum(as.numeric(sim.stat[i, 
                  ] >= obs.stat[i])))/(NREP + 1)
            else pvalue[i] <- (1 + sum(as.numeric(sim.stat[i] >= 
                obs.stat[i])))/(NREP + 1)
       }
      summary <- matrix(c(lags,ans[,2],ans[,3],pvalue),ncol=4)
    dimnames(summary) <- list(rep("", length(lags)),c("Lags","Statistic","df","p-value"))
  return(summary)
    }
}
