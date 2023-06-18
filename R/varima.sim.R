"varima.sim" <-
function (model=list(ar=NULL,ma=NULL,d=NULL,sar=NULL,sma=NULL,D=NULL,period=NULL),
n, k = 1, constant = NA, trend = NA, demean = NA, innov = NULL,
innov.dist=c("Gaussian","t","bootstrap"),...)
{
    	dots <- list(...)
    	trunc.lag <- dots$trunc.lag            ## truncation infinite series
    	dft <- dots$dft                        ## df for t-distribution
    	sigma <- dots$sigma
    	innov.dist <- match.arg(innov.dist)
      ar <- model$ar
      ma <- model$ma
      d <- model$d
      ar.season <- model$sar
      ma.season <- model$sma
      D <- model$D
      period <- model$period
	    if (!is.null(ar) && !inherits(ar,"array") && !inherits(ar,"numeric"))
		stop("ar must be enterd as NULL or numeric or array with dimension (k*k*p)")
	    if (!is.null(ma) && !inherits(ma,"array") && !inherits(ma,"numeric"))
		stop("ma must be enterd as NULL or numeric or array with dimension (k*k*q)")
	    if (!is.null(ar.season) && !inherits(ar.season,"array") && !inherits(ar.season,"numeric"))
		stop("ar.season must be enterd as NULL or numeric or array with dimension (k*k*ps)")
	    if (!is.null(ma.season) && !inherits(ma.season,"array") && !inherits(ma.season,"numeric"))
		stop("ma.season must be enterd as NULL or numeric or array with dimension (k*k*qs)")
          if (is.null(sigma)) sigma <- diag(k)
	    else sigma <- as.matrix(sigma,nrow=k)
	    if (all(is.null(d))||all(d==0))
	        d <- rep(0, k)
	    if (length(d) != k)
	      stop("d must be entered as a vector with length equal to number of sigma rows")
 	    if (any(d < 0))
            stop("number of differences must be a nonnegative integer/integers")
	    if (all(is.null(D))||all(D==0))
	        D<- rep(0, k)
	    if (length(D) != k)
	      stop("D must be entered as a vector with length equal to number of sigma rows")
 	    if (any(D < 0))
            stop("number of season differences must be a nonnegative integer/integers")
	    if (all(is.na(constant)))
 	        constant <- rep(0, k)
	    if (length(constant) != k)
	      stop("constant must be entered as a vector with length equal to number of sigma rows")
	    if (all(is.na(trend)))
      	  trend <- rep(0, k)
	    if (length(trend) != k)
	        stop("trend must be entered as a vector with length equal to number of sigma rows")
	    if (all(is.na(demean)))
	        demean <- rep(0, k)
	    if (length(demean) != k)
	        stop("demean must be entered as a vector with length equal to number of sigma rows")
	    if (inherits(ar,"numeric"))
	        ar <- array(ar,dim=c(1,1,length(ar)))
	    if (inherits(ma,"numeric"))
	        ma <- array(ma,dim=c(1,1,length(ma)))
	    if (inherits(ar.season,"numeric"))
	        ar.season <- array(ar.season,dim=c(1,1,length(ar.season)))
	    if (inherits(ma.season,"numeric"))
	        ma.season <- array(ma.season,dim=c(1,1,length(ma.season)))
 	    if (all(ar == 0))
 	        ar <- NULL
	    if (all(ma == 0))
              ma <- NULL
 	    if (all(ar.season == 0))
 	        ar.season <- NULL
	    if (all(ma.season == 0))
              ma.season <- NULL
          p <- ifelse(is.null(ar),0,dim(ar)[3])
          q <- ifelse(is.null(ma),0,dim(ma)[3])
          ps <- ifelse(is.null(ar.season),0,dim(ar.season)[3])
          qs <- ifelse(is.null(ma.season),0,dim(ma.season)[3])
          if (p > 0 && ((dim(ar)[1] != dim(ar)[2])||dim(ar)[2] != k))
             stop("Incompatible dimensions of the ar coefficients and k")
          if (q > 0 && ((dim(ma)[1] != dim(ma)[2])||dim(ma)[2] != k))
             stop("Incompatible dimensions of the ma coefficients and k")
          if (ps > 0 && ((dim(ar.season)[1] != dim(ar.season)[2])||dim(ar.season)[2] != k))
             stop("Incompatible dimensions of the sar coefficients and k")
          if (qs > 0 && ((dim(ma.season)[1] != dim(ma.season)[2])||dim(ma.season)[2] != k))
             stop("Incompatible dimensions of the sma coefficients and k")
          if ((ps > 0 && is.null(period)) || (qs && is.null(period)))
             period <- 12
          if(!is.null(innov)){ ## initial innovations (signal processing) is given
                innov <- as.matrix(innov)
                 stopifnot (NROW(innov) == n && NCOL(innov) == k)
          }
       s <- period
       if (ps != 0){  ## Big.ar = ar(B) * ar(B^s)
           ar.Season <- array(c(diag(k),-ar.season),dim=c(k,k,ps+1))
           Extend.ar.season <- array(numeric(k*k*(k+ps*s)),dim=c(k,k,1+ps*s))
          for (i in (1:(ps + 1)))
            Extend.ar.season[,,1+s*{i-1}] <- ar.Season[,,i]
          if (p !=0){ ## ar = ar(B) * ar(B^s)
         	 Big.ar <- array(numeric(k*k*(p+1)*(k+p+ps*s)),dim=c(k,k*(p+1),1+p+ps*s))
         	 Big.ar[,((k*p+1):((p+1)*k)),1:(dim(Extend.ar.season)[3])] <- Extend.ar.season
        	   for(i in 1:p){
        	      for(j in 1:(dim(Extend.ar.season)[3]))
        	     Big.ar[,((1+(i-1)*k):(i*k)),i+j] <- crossprod(t(-ar[,,i]), Extend.ar.season[,,j])
       	   }
             out.Big.ar <- array(numeric(k*k*(1+p+ps*s)),dim=c(k,k,1+p+ps*s))
	       for (i in 1:(1+p+ps*s)){
               Mat1 <- matrix(rep(0,k*k),k,k)
                for (j in (1:(p+1)))
                  Mat1 <- Big.ar[,((1+(j-1)*k):(j*k)),i] + Mat1
               out.Big.ar[,,i] <- Mat1
              }
           ar <- -array(out.Big.ar[,,-1],dim=c(k,k,p+ps*s))
          }
          else { ## pure seasonal series: ar = ar(B^s)
             ar <- -array(Extend.ar.season[,,-1],dim=c(k,k,ps*s))
          }
         p <- p + ps*s
        }
        else if (qs != 0) { ## Big.ma = ma(B) * ma(B^s)
           ma.Season <- array(c(diag(k),-ma.season),dim=c(k,k,qs+1))
           Extend.ma.season <- array(numeric(k*k*(k+qs*s)),dim=c(k,k,1+qs*s))
          for (i in (1:(qs + 1)))
            Extend.ma.season[,,1+s*{i-1}] <- ma.Season[,,i]
          if (q !=0){ ## ma = ma(B) * ma(B^s)
         	 Big.ma <- array(numeric(k*k*(q+1)*(k+q+qs*s)),dim=c(k,k*(q+1),1+q+qs*s))
         	 Big.ma[,((k*q+1):((q+1)*k)),1:(dim(Extend.ma.season)[3])] <- Extend.ma.season
        	   for(i in 1:q){
        	      for(j in 1:(dim(Extend.ma.season)[3]))
        	     Big.ma[,((1+(i-1)*k):(i*k)),i+j] <- crossprod(t(-ma[,,i]), Extend.ma.season[,,j])
       	   }
             out.Big.ma <- array(numeric(k*k*(1+q+qs*s)),dim=c(k,k,1+q+qs*s))
	       for (i in 1:(1+q+qs*s)){
               Mat2 <- matrix(rep(0,k*k),k,k)
                for (j in (1:(q+1)))
                  Mat2 <- Big.ma[,((1+(j-1)*k):(j*k)),i] + Mat2
               out.Big.ma[,,i] <- Mat2
              }
           ma <- -array(out.Big.ma[,,-1],dim=c(k,k,q+qs*s))
          }
          else { ## pure seasonal series: ma = ma(B^s)
             ma <- -array(Extend.ma.season[,,-1],dim=c(k,k,qs*s))
          }
        q <- q + qs*s
     }
     if (p == 0) { ## autoregressive part is not included in VARMA(p,q)
            if(innov.dist == "Gaussian") ## stimulate innovation from normal distribution
                 epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(n + q)),ncol=n + q)))
            else if(innov.dist == "t"){ ## stimulate innovation from t-distribution
              if (is.null(dft)|| dft <= 0)
                stop("degrees of freedom for t-distribution is missing or entered as negative value")
              else
                epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(n + q)),ncol=n + q)))/sqrt(rchisq(k*(n + q),dft)/dft)
             }
            else if (innov.dist == "bootstrap" && is.null(innov))
                stop("innovation series needed for bootstraping error term is missing!")
            else if (innov.dist == "bootstrap" && !is.null(innov)) ## stimulate innovation using resampled errors instead of distribution
                epsilon <- ts(matrix(sample(x=innov, size=k*(n+q), replace=TRUE),ncol = k))
         if (q == 0){  ## Simulate white noise - VARMA(0,0)
                 Sloptrend <- t(matrix(trend * rep(1:NROW(epsilon),each=k),nrow=k,ncol=NROW(epsilon)))
                 trendTerm <- sweep(Sloptrend, 2L, -constant, check.margin = FALSE)
                 CenterData <- scale(epsilon, center = -demean, scale = FALSE)
                 Sim.Series <- trendTerm + CenterData
                 for (i in 1:length(d)){
                   if(d[i]>0)
                     Sim.Series[,i] <- as.matrix(diffinv(Sim.Series[,i], differences = d[i]))[-(1:d[i]),]
                   else
                     Sim.Series[,i] <- Sim.Series[,i]
                   }
                 for (j in 1:length(D)){
                   if(D[j]>0)
                     Sim.Series[,j] <- as.matrix(diffinv(Sim.Series[,j], lag = s, differences = D[j]))[-(1:(s*D[j])),]
                   else
                     Sim.Series[,j] <- Sim.Series[,j]
                   }
                 return(ts(Sim.Series))
           }
           else ## Simulate VMA(q)
                 InvertQ(ma)
                 psi <- array(c(diag(k), -ma), dim = c(k, k, q + 1))
            if(innov.dist == "Gaussian") ## stimulate innovation from normal distribution
                 epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(n + q)),ncol=n + q)))
            else if(innov.dist == "t"){ ## stimulate innovation from t-distribution
              if (is.null(dft)|| dft <= 0)
                stop("degrees of freedom for t-distribution is missing or entered as negative value")
              else
                epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(n + q)),ncol=n + q)))/sqrt(rchisq(k*(n + q),dft)/dft)
             }
           else if (innov.dist == "bootstrap" && is.null(innov))
                stop("innovation series needed for bootstraping error term is missing!")
            else if (innov.dist == "bootstrap" && !is.null(innov)) ## stimulate innovation using resampled errors instead of distribution
                epsilon <- ts(matrix(sample(x=innov, size=k*(n+q), replace=TRUE),ncol = k))
                 Sim.VMA <- vma.sim(psi = psi, a = epsilon)
                 Sloptrend <- t(matrix(trend * rep(1:NROW(Sim.VMA),each=k),nrow=k,ncol=NROW(Sim.VMA)))
                 trendTerm <- sweep(Sloptrend, 2L, -constant, check.margin = FALSE)
                 CenterData <- scale(Sim.VMA, center = -demean, scale = FALSE)
                 Sim.Series <- trendTerm + CenterData
                 for (i in 1:length(d)){
                    if(d[i]>0)
                      Sim.Series[,i] <- as.matrix(diffinv(Sim.Series[,i], differences = d[i]))[-(1:d[i]),]
                    else
                      Sim.Series[,i] <- Sim.Series[,i]
                 }
                 for (j in 1:length(D)){
                   if(D[j]>0)
                     Sim.Series[,j] <- as.matrix(diffinv(Sim.Series[,j], lag = s, differences = D[j]))[-(1:(s*D[j])),]
                   else
                     Sim.Series[,j] <- Sim.Series[,j]
                 }
                return(ts(Sim.Series))
              }
          else ## autoregressive part is included in VARMA(p,q)
            if (is.null(trunc.lag))
               trunc.lag <- min(100,ceiling(n/3))
            FirstSim.Series <- matrix(numeric(0), nrow = n, ncol = k)
            r <- max(p, q)
            psi <- ImpulseVMA(ar = ar, ma = ma, trunc.lag = trunc.lag)
            if(innov.dist == "Gaussian")
                epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(trunc.lag + r)),ncol=trunc.lag + r)))
            else if(innov.dist == "t"){ ## stimulate innovation from t-distribution
              if (is.null(dft)|| dft <= 0)
                stop("degrees of freedom for t-distribution is missing or entered as negative value")
              else
                epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(trunc.lag + r)),ncol=trunc.lag + r)))/sqrt(rchisq(k*(trunc.lag + r),dft)/dft)
              }
            else if (innov.dist == "bootstrap" && is.null(innov))
                stop("innovation series needed for bootstraping error term is missing!")
            else if (innov.dist == "bootstrap" && !is.null(innov))
                epsilon <- matrix(sample(x=innov, size=k*(trunc.lag + r), replace=TRUE),ncol = k)
            FirstSim.Series[1:r, ] <- vma.sim(psi = psi, a = epsilon)
            a <- matrix(epsilon[1:r, ], nrow = r, ncol = k)
            if(innov.dist == "Gaussian")
               epsilon <- rbind(a,t(crossprod(chol(sigma),matrix(stats::rnorm(k*n),ncol=n))))
            else if(innov.dist == "t"){
              if (is.null(dft)|| dft <= 0)
                stop("degrees of freedom for t-distribution is missing or entered as negative value")
               epsilon <- rbind(a,t(crossprod(chol(sigma),matrix(stats::rnorm(k*n),ncol=n)))/sqrt(rchisq(k*n,dft)/dft))
            }
             else if (innov.dist == "bootstrap" && is.null(innov))
                stop("innovation series needed for bootstraping error term is missing!")
            else if (innov.dist == "bootstrap" && !is.null(innov))
                epsilon <- rbind(a,matrix(sample(x=innov,size=k*n,replace=TRUE),ncol=k))
            if (q > 0) { ## Simulate VARMA(p,q)
                extend.psi <- array(c(diag(k), -ma, rep(0, k * k *(n - q))), dim = c(k, k, n + 1))
                u <- matrix(numeric(0), nrow = n, ncol = k)
                for (i in (q + 1):(n + q)) {
                   out <- 0
                   for (j in 0:q) {
                      out = out + crossprod(t(extend.psi[, , j + 1]),epsilon[i - j, ])
                   }
                  u[i - q, ] <- out
                }
             }
            else
               u <- epsilon  ## Simulate VAR(p)
               for (i in (r + 1):n) {
                 temp2 <- 0
                for (j in 1:p)
                     temp2 <- temp2 + crossprod(t(ar[, , j]),FirstSim.Series[i - j, ])
                FirstSim.Series[i, ] <- temp2 + u[i, ]
               }
              Sloptrend <- t(matrix(trend * rep(1:NROW(FirstSim.Series),each=k),nrow=k,ncol=NROW(FirstSim.Series)))
              trendTerm <- sweep(Sloptrend, 2L, -constant, check.margin = FALSE)
              CenterData <- scale(FirstSim.Series, center = -demean, scale = FALSE)
              Sim.Series <- trendTerm + CenterData
            for (i in 1:length(d)){
              if(d[i]>0)
                 Sim.Series[,i] <- as.matrix(diffinv(Sim.Series[,i], differences = d[i]))[-(1:d[i]),]
              else
                 Sim.Series[,i] <- Sim.Series[,i]
             }
            for (j in 1:length(D)){
              if(D[j]>0)
                 Sim.Series[,j] <- as.matrix(diffinv(Sim.Series[,j], lag = s, differences = D[j]))[-(1:(s*D[j])),]
              else
                 Sim.Series[,j] <- Sim.Series[,j]
            }
          return(ts(Sim.Series))
}
