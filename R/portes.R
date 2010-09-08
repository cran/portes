"portes" <-
  function(obj, lags=seq(5,30,5), statistic=c("GVStat","LiMcLeod"), order=0, 
   SpawnSlaves=c("default","Rmpi","snow"), nslaves=2, NREP=1000, InfiniteVarianceQ = FALSE,
       SquaredQ=FALSE,Set.Seed=TRUE)
{
     statistic <- match.arg(statistic)
     SpawnSlaves <- match.arg(SpawnSlaves)
           TestType<-"0"
        if (all(class(obj)=="ts")||all(class(obj)=="numeric")||all(class(obj)=="matrix")||all((class(obj)[1]=="mts" && class(obj)[2]=="ts"))) 
           TestType <- "1"  #bootstrap test for randomness 
        else if (all(class(obj)=="ar")||all(class(obj)=="Arima")||all(class(obj)=="arima0")||all(class(obj)=="FitAR") || all(class(obj)=="FitFGN"))
           TestType <- "2"  #bootstrap VAR or arima diagnostic check
        if (TestType == "0")
           stop("obj must be class ar, Arima, arima0, FitAR, FitFGN, ts, numeric, matrix, or (mts ts)")
        if (statistic=="GVStat") 
    		obs.stat <- GVStat(obj,lags,order,SquaredQ)[,2]
        else 
            obs.stat <- LiMcLeod(obj,lags,order,SquaredQ)[,2]
                    if (all(class(obj) == "ar")){ 
	               p <- obj$order
                     q <- 0
                     Order <- p+q
                     res <- ts(as.matrix(obj$resid)[-(1:Order),])
	               n <- NROW(res)
	               k <- NCOL(res)
                       if (k==1) 
                          Model <- 1
                       else
                          Model <- 2
                        sigma <- obj$var.pred
                            if (is.array(obj$ar)){
                               arrayphi <- array(numeric(k * k * p), dim = c(k^2, p))
                              for (i in 1:p) arrayphi[, i] <- c(obj$ar[i, , ])
                               phi <- array(c(arrayphi), dim = c(k, k, p))
                            }
                            else
                               phi <- obj$ar
                            theta <- NULL
                           if (!is.null(obj$x.intercept)) 
                                intercept <- obj$x.intercept
                           else 
                                intercept <- rep(0,k)
                     }
                     else if (all(class(obj) == "Arima")||all(class(obj) == "arima0")){
                         Model <- 3
                         res <- obj$resid
                         n <- length(res)
                         k <- 1
                         pq <- eval(parse(text = (as.character(obj$call))[[3]]))
                         p <- pq[1]
                         q <- pq[3]
                         Order <- p+q
                         phi <- theta <- NULL
                         if (p > 0) 
                             phi <- obj$model$phi
                         if (q > 0) 
                             theta <- obj$model$theta
                          sigma <- obj$sigma2
                          intercept <- 0
                     }
                     else if (all(class(obj) == "FitAR")){
                          Model <- 4
                          res <- obj$res
                          n <- length(res)
                          k <- 1
                          p <- length(obj$phiHat)
                          q <- 0
                          Order <- p+q
                          phi <- obj$phiHat
                          theta <- NULL
                          sigma <- obj$sigsqHat
                          intercept <- obj$muHat
                     }
                     else if (all(class(obj) == "FitFGN")){
                           Model <- 5
                           phi <- theta <- NULL
                           res <- obj$res
                           n <- length(res)
                           k <- 1
                           p <- 1
                           q <- 0
                           Order <- 1
                           n <- length(obj$res)
                           H <- obj$H
                           sigma <- sqrt(obj$sigsq)
                           intercept <- obj$muHat
                     }
                     else {
                         Order <- order
                         res <- ts(as.matrix(obj)) 
		             n <- NROW(res)
                         k <- NCOL(res)
    	                   sigma <- matrix(acf(res, lag.max = 1, plot = FALSE, type = "covariance")$acf[1,,],k,k)
                     }
             if (InfiniteVarianceQ){
                 StableParameters<-matrix(FitStable(res),ncol=4)
                     ALPHA <- StableParameters[,1] #shape
                     BETA <- StableParameters[,2] #skewness
                     GAMMA <- StableParameters[,3] #scale
                     DELTA <- StableParameters[,4] #location
	          ALPHA <<- ALPHA
	          BETA <<- BETA
	          GAMMA <<- GAMMA
	          DELTA <<- DELTA
              }
             else StableParameters <- NA
             Trunc.Series <- min(100,ceiling(n/3))
          if (TestType=="1"){
	         (OneBoot <- function(){
                        if (InfiniteVarianceQ)
                           rboot<- ts(rstable(n, ALPHA, BETA, GAMMA, DELTA))
                        else
                           rboot<- t(crossprod(chol(sigma),matrix(rnorm(k*n),ncol=n)))
                        if(statistic=="GVStat")
                           OneSim.stat<- GVStat(rboot,lags,Order,SquaredQ)[,2]
                        else 
                           OneSim.stat<- LiMcLeod(rboot,lags,Order,SquaredQ)[,2]
                       return(OneSim.stat)
                  }
                )
           }
           else{
		     (OneBoot <- function(){
  			    if (Model==1){
     				     if (InfiniteVarianceQ)
      				  Sim.Data <- arima.sim(n=n,list(ar=phi, ma=theta),rand.gen=function(n, ...) ts(intercept+rstable(n, ALPHA, BETA, GAMMA, DELTA)))
                              else
                                Sim.Data <- arima.sim(n = n, list(ar=phi, ma=theta),sd = sqrt(sigma),mean=intercept)
                             FitSimModel <- ar(Sim.Data, aic = FALSE, order.max = p)
                             rboot <- FitSimModel$resid[-(1:p)]
                             if(statistic=="GVStat")
                                  OneSim.stat<- GVStat(rboot,lags,Order,SquaredQ)[,2]
                             else 
                                  OneSim.stat<- LiMcLeod(rboot,lags,Order,SquaredQ)[,2]
                            return(OneSim.stat)
                       }
                       else if (Model==2){
                                Sim.Data <- simvarma(phi = phi, theta = theta, sigma = sigma, intercept = intercept, n = n, StableParameters = StableParameters, Trunc.Series = Trunc.Series)
                           if(!all(intercept==0))                           
                             FitSimModel <- ar.ols(Sim.Data, aic = FALSE, intercept = TRUE, order.max = p)
                           else
                             FitSimModel <- ar.ols(Sim.Data, aic = FALSE, intercept = FALSE, order.max = p)
                             rboot <- ts(as.matrix(FitSimModel$resid)[-(1:p), ])
                             if(statistic=="GVStat")
                                  OneSim.stat<- GVStat(rboot,lags,Order,SquaredQ)[,2]
                             else 
                                  OneSim.stat<- LiMcLeod(rboot,lags,Order,SquaredQ)[,2]
                           return(OneSim.stat)
                       }
                       else if (Model==3){
                               if (InfiniteVarianceQ)
                                  Sim.Data <- arima.sim(n=n,list(ar=phi, ma=theta),rand.gen=function(n, ...) ts(intercept+rstable(n, ALPHA, BETA, GAMMA, DELTA)))
                                else
                                  Sim.Data <- arima.sim(n = n, list(ar=phi, ma=theta),sd = sqrt(sigma),mean=intercept)
                               FitSimModel <- arima(Sim.Data, c(p,0,q))
                               rboot <- FitSimModel$resid
                             if(statistic=="GVStat")
                                  OneSim.stat<- GVStat(rboot,lags,Order,SquaredQ)[,2]
                             else 
                                  OneSim.stat<- LiMcLeod(rboot,lags,Order,SquaredQ)[,2]
                            return(OneSim.stat)
                         }
                         else if (Model==4){
                               if (InfiniteVarianceQ)
                                  Sim.Data <- arima.sim(n=n,list(ar=phi, ma=theta),rand.gen=function(n, ...) ts(intercept+rstable(n, ALPHA, BETA, GAMMA, DELTA)))
                                else
                                  Sim.Data <- arima.sim(n = n, list(ar=phi, ma=theta),sd = sqrt(sigma),mean=intercept)
                               FitSimModel <- FitAR(Sim.Data, p)
                               rboot <- ts(FitSimModel$res)
                             if(statistic=="GVStat")
                                  OneSim.stat<- GVStat(rboot,lags,Order,SquaredQ)[,2]
                             else 
                                  OneSim.stat<- LiMcLeod(rboot,lags,Order,SquaredQ)[,2]
                            return(OneSim.stat)
                          }
                          else if (Model==5){
                                 Sim.Data <- intercept + SimulateFGN(n, H) * sigma
                              FitSimModel <- FitFGN(Sim.Data)
                              rboot <- ts(FitSimModel$res)
                             if(statistic=="GVStat")
                                  OneSim.stat<- GVStat(rboot,lags,Order,SquaredQ)[,2]
                             else 
                                  OneSim.stat<- LiMcLeod(rboot,lags,Order,SquaredQ)[,2]
                            return(OneSim.stat)
                          }
                      }
                  )
           }
         if (SpawnSlaves=="default"){ 
             if (Set.Seed) 
                 set.seed(21597341) 
             sim.stat <- replicate(NREP,OneBoot())
          }
         else if (SpawnSlaves=="Rmpi"){
              mpi.spawn.Rslaves(nslaves = nslaves)
                 if (Set.Seed) 
                  mpi.setup.rngstream(21597341)
            mpi.bcast.cmd(library(akima))
            mpi.bcast.cmd(library(FitAR))
            mpi.bcast.cmd(library(FGN)) 
	      mpi.bcast.Robj2slave(k)
	      mpi.bcast.Robj2slave(n)
	      mpi.bcast.Robj2slave(Order)
	      mpi.bcast.Robj2slave(obj)
	      mpi.bcast.Robj2slave(lags)
	      mpi.bcast.Robj2slave(NREP)
	      mpi.bcast.Robj2slave(res)
	      mpi.bcast.Robj2slave(sigma)
	      mpi.bcast.Robj2slave(InfiniteVarianceQ)
	      mpi.bcast.Robj2slave(SquaredQ)
	      mpi.bcast.Robj2slave(Trunc.Series)
              if (InfiniteVarianceQ){
	         mpi.bcast.Robj2slave(ALPHA)
	         mpi.bcast.Robj2slave(BETA)
	         mpi.bcast.Robj2slave(GAMMA)
	         mpi.bcast.Robj2slave(DELTA)
	         mpi.bcast.Robj2slave(interpp)
	         mpi.bcast.Robj2slave(interpp.old)
	         mpi.bcast.Robj2slave(FitStable)
	         mpi.bcast.Robj2slave(rstable)
               }
              if (TestType=="2"){
	         mpi.bcast.Robj2slave(p)
	         mpi.bcast.Robj2slave(Model)
	         mpi.bcast.Robj2slave(phi)
	         mpi.bcast.Robj2slave(theta)
               mpi.bcast.Robj2slave(intercept)
              }
	       mpi.bcast.Robj2slave(statistic)
	       mpi.bcast.Robj2slave(blockToeplitz)
	       mpi.bcast.Robj2slave(Get.Resid)
	       mpi.bcast.Robj2slave(GVStat)
	       mpi.bcast.Robj2slave(ImpulseVMA)
	       mpi.bcast.Robj2slave(InvertQ)
	       mpi.bcast.Robj2slave(LiMcLeod)
	       mpi.bcast.Robj2slave(OneBoot)
	       mpi.bcast.Robj2slave(simvarma)
	       mpi.bcast.Robj2slave(simvma)
           sim.stat <- mpi.parReplicate(NREP,OneBoot())
          mpi.close.Rslaves()
         }
         else if (SpawnSlaves=="snow"){ 
	      n <<- n
	      k <<- k
            Order <<- Order
	      obj <<- obj
	      lags <<- lags
	      res <<- res
            NREP <<- NREP
            if (TestType=="2"){
	      p <<- p
	      Model <<- Model
            phi <<- phi
            theta <<- theta
            intercept <<- intercept
           }
            sigma <<- sigma
	      SquaredQ <<- SquaredQ
	      statistic <<- statistic
	      InfiniteVarianceQ <<- InfiniteVarianceQ
	      Trunc.Series <<- Trunc.Series
	      Set.Seed <<- Set.Seed
            OneBoot <<- OneBoot 
            cl <- makeCluster(nslaves)
               if (Set.Seed) 
                   clusterSetupRNG(cl,seed=21597341)
               clusterEvalQ(cl,library(FitAR))
               clusterEvalQ(cl,library(FGN))
               clusterEvalQ(cl,library(akima))
               clusterExport(cl, list("k", "n","obj", "Order","sigma","lags","NREP","res","InfiniteVarianceQ","SquaredQ","Trunc.Series"))
               clusterExport(cl, list("blockToeplitz","Get.Resid","GVStat","ImpulseVMA","InvertQ","LiMcLeod","OneBoot","simvarma","simvma","statistic"))
             if (InfiniteVarianceQ)
                 clusterExport(cl, list("ALPHA", "BETA", "GAMMA", "DELTA", "interpp", "interpp.old","MCTable3","MCTable4","MCTable5","MCTable7","FitStable","rstable"))    
             if (TestType=="2")
                 clusterExport(cl, list("p","Model","phi","theta","intercept"))
             sim.stat <- parSapply(cl, 1:NREP, function(j)(OneBoot()))
            stopCluster(cl)
           }
           pvalue <- numeric(length(lags))
              for (i in 1:length(lags))
                  pvalue[i] <- (1+sum(as.numeric(sim.stat[i,]>=obs.stat[i])))/(NREP+1)
              summary <- matrix(c(lags,pvalue),ncol=2)
              dimnames(summary) <- list(rep("", length(lags)),c("Lags","p-value"))
          return(summary)
}
