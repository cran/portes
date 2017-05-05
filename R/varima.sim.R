"varima.sim" <-
function (model=list(ar=NULL,ma=NULL,d=NULL,sar=NULL,sma=NULL,D=NULL,period=NULL),
n, k = 1, constant = NA, trend = NA, demean = NA, innov = NULL, 
innov.dist=c("Gaussian","t","stable","bootstrap"),...) 
{
    	dots <- list(...) 
    	trunc.lag <- dots$trunc.lag            ## truncation infinite series 
    	dft <- dots$dft                        ## df for t-distribution
    	par.stable <- dots$par.stable          ## stable 4 parameters 
    	sigma <- dots$sigma
    	innov.dist <- match.arg(innov.dist)
      phi <- model$ar
      theta <- model$ma
      d <- model$d
      phi.season <- model$sar
      theta.season <- model$sma
      D <- model$D
      period <- model$period
	    if (!is.null(phi) && class(phi)!="array" && class(phi)!="numeric")
		stop("phi must be enterd as NULL or numeric or array with dimension (k*k*p)")
	    if (!is.null(theta) && class(theta)!="array" && class(theta)!="numeric")
		stop("theta must be enterd as NULL or numeric or array with dimension (k*k*q)")
	    if (!is.null(phi.season) && class(phi.season)!="array" && class(phi.season)!="numeric")
		stop("phi.season must be enterd as NULL or numeric or array with dimension (k*k*ps)")
	    if (!is.null(theta.season) && class(theta.season)!="array" && class(theta.season)!="numeric")
		stop("theta.season must be enterd as NULL or numeric or array with dimension (k*k*qs)")
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
	    if (class(phi) == "numeric")
	        phi <- array(phi,dim=c(1,1,length(phi)))
	    if (class(theta)=="numeric")
	        theta <- array(theta,dim=c(1,1,length(theta)))
	    if (class(phi.season) == "numeric")
	        phi.season <- array(phi.season,dim=c(1,1,length(phi.season)))
	    if (class(theta.season)=="numeric")
	        theta.season <- array(theta.season,dim=c(1,1,length(theta.season)))
 	    if (all(phi == 0))
 	        phi <- NULL
	    if (all(theta == 0))
              theta <- NULL
 	    if (all(phi.season == 0))
 	        phi.season <- NULL
	    if (all(theta.season == 0))
              theta.season <- NULL
          p <- ifelse(is.null(phi),0,dim(phi)[3])
          q <- ifelse(is.null(theta),0,dim(theta)[3])
          ps <- ifelse(is.null(phi.season),0,dim(phi.season)[3])
          qs <- ifelse(is.null(theta.season),0,dim(theta.season)[3])
          if (p > 0 && ((dim(phi)[1] != dim(phi)[2])||dim(phi)[2] != k))
             stop("Incompatible dimensions of the ar coefficients and k")
          if (q > 0 && ((dim(theta)[1] != dim(theta)[2])||dim(theta)[2] != k))
             stop("Incompatible dimensions of the ma coefficients and k")
          if (ps > 0 && ((dim(phi.season)[1] != dim(phi.season)[2])||dim(phi.season)[2] != k))
             stop("Incompatible dimensions of the sar coefficients and k")
          if (qs > 0 && ((dim(theta.season)[1] != dim(theta.season)[2])||dim(theta.season)[2] != k))
             stop("Incompatible dimensions of the sma coefficients and k")
          if ((ps > 0 && is.null(period)) || (qs && is.null(period)))
             period <- 12
          if(!is.null(innov)){ ## initial innovations (signal processing) is given
                innov <- as.matrix(innov)
                 stopifnot (NROW(innov) == n && NCOL(innov) == k)
          }
          if(!is.null(par.stable)) {
            par.stable <- matrix(par.stable,ncol=4)
              if (dim(par.stable)[1]!=k)
                  stop("par.stable must be a numeric vector/matrix of dimension k-by-4 with 4 stable parameters: Alpha, Beta, Scale, and Location")
              Alpha<-par.stable[,1]
              Beta<-par.stable[,2]
              Scale<-par.stable[,3]
              Location<-par.stable[,4]
          }
       s <- period
       if (ps != 0){  ## Big.Phi = phi(B) * Phi(B^s)
           Phi.Season <- array(c(diag(k),-phi.season),dim=c(k,k,ps+1))
           Extend.phi.season <- array(numeric(k*k*(k+ps*s)),dim=c(k,k,1+ps*s))
          for (i in (1:(ps + 1)))
            Extend.phi.season[,,1+s*{i-1}] <- Phi.Season[,,i]
          if (p !=0){ ## phi = phi(B) * Phi(B^s)
         	 Big.Phi <- array(numeric(k*k*(p+1)*(k+p+ps*s)),dim=c(k,k*(p+1),1+p+ps*s))
         	 Big.Phi[,((k*p+1):((p+1)*k)),1:(dim(Extend.phi.season)[3])] <- Extend.phi.season
        	   for(i in 1:p){
        	      for(j in 1:(dim(Extend.phi.season)[3]))
        	     Big.Phi[,((1+(i-1)*k):(i*k)),i+j] <- crossprod(t(-phi[,,i]), Extend.phi.season[,,j])
       	   }
             out.Big.Phi <- array(numeric(k*k*(1+p+ps*s)),dim=c(k,k,1+p+ps*s))
	       for (i in 1:(1+p+ps*s)){
               Mat1 <- matrix(rep(0,k*k),k,k)
                for (j in (1:(p+1))) 
                  Mat1 <- Big.Phi[,((1+(j-1)*k):(j*k)),i] + Mat1   
               out.Big.Phi[,,i] <- Mat1
              }
           phi <- -array(out.Big.Phi[,,-1],dim=c(k,k,p+ps*s))
          }
          else { ## pure seasonal series: phi = Phi(B^s)
             phi <- -array(Extend.phi.season[,,-1],dim=c(k,k,ps*s))
          }
         p <- p + ps*s
        }
        else if (qs != 0) { ## Big.Theta = theta(B) * Theta(B^s)
           Theta.Season <- array(c(diag(k),-theta.season),dim=c(k,k,qs+1))
           Extend.theta.season <- array(numeric(k*k*(k+qs*s)),dim=c(k,k,1+qs*s))
          for (i in (1:(qs + 1)))
            Extend.theta.season[,,1+s*{i-1}] <- Theta.Season[,,i]
          if (q !=0){ ## theta = theta(B) * Theta(B^s)
         	 Big.Theta <- array(numeric(k*k*(q+1)*(k+q+qs*s)),dim=c(k,k*(q+1),1+q+qs*s))
         	 Big.Theta[,((k*q+1):((q+1)*k)),1:(dim(Extend.theta.season)[3])] <- Extend.theta.season
        	   for(i in 1:q){
        	      for(j in 1:(dim(Extend.theta.season)[3]))
        	     Big.Theta[,((1+(i-1)*k):(i*k)),i+j] <- crossprod(t(-theta[,,i]), Extend.theta.season[,,j])
       	   }
             out.Big.Theta <- array(numeric(k*k*(1+q+qs*s)),dim=c(k,k,1+q+qs*s))
	       for (i in 1:(1+q+qs*s)){
               Mat2 <- matrix(rep(0,k*k),k,k)
                for (j in (1:(q+1))) 
                  Mat2 <- Big.Theta[,((1+(j-1)*k):(j*k)),i] + Mat2   
               out.Big.Theta[,,i] <- Mat2
              }
           theta <- -array(out.Big.Theta[,,-1],dim=c(k,k,q+qs*s))
          }
          else { ## pure seasonal series: theta = Theta(B^s)
             theta <- -array(Extend.theta.season[,,-1],dim=c(k,k,qs*s))
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
            else if(innov.dist == "stable") ## stimulate innovation from stable distribution         
                epsilon <- rStable(n=n + q, Alpha=Alpha,Beta=Beta,Scale=Scale,Location=Location)
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
                 InvertQ(theta) 
                 psi <- array(c(diag(k), -theta), dim = c(k, k, q + 1))
            if(innov.dist == "Gaussian") ## stimulate innovation from normal distribution
                 epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(n + q)),ncol=n + q)))
            else if(innov.dist == "t"){ ## stimulate innovation from t-distribution
              if (is.null(dft)|| dft <= 0) 
                stop("degrees of freedom for t-distribution is missing or entered as negative value") 
              else
                epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(n + q)),ncol=n + q)))/sqrt(rchisq(k*(n + q),dft)/dft)
             }
            else if(innov.dist == "stable")  ## stimulate innovation from stable distribution            
                epsilon <- rStable(n=n + q, Alpha=Alpha,Beta=Beta,Scale=Scale,Location=Location)
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
            psi <- ImpulseVMA(phi = phi, theta = theta, trunc.lag = trunc.lag)
            if(innov.dist == "Gaussian")
                epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(trunc.lag + r)),ncol=trunc.lag + r)))
            else if(innov.dist == "t"){ ## stimulate innovation from t-distribution
              if (is.null(dft)|| dft <= 0) 
                stop("degrees of freedom for t-distribution is missing or entered as negative value") 
              else
                epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(trunc.lag + r)),ncol=trunc.lag + r)))/sqrt(rchisq(k*(trunc.lag + r),dft)/dft)
              }            
            else if(innov.dist == "stable") ## stimulate innovation from stable distribution                
                epsilon <- rStable(n=trunc.lag + r, Alpha=Alpha,Beta=Beta,Scale=Scale,Location=Location)
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
            else if(innov.dist == "stable")   
               epsilon <- rbind(a,rStable(n=n,Alpha=Alpha,Beta=Beta,Scale=Scale,Location=Location))
            else if (innov.dist == "bootstrap" && is.null(innov)) 
                stop("innovation series needed for bootstraping error term is missing!")
            else if (innov.dist == "bootstrap" && !is.null(innov))
                epsilon <- rbind(a,matrix(sample(x=innov,size=k*n,replace=TRUE),ncol=k))
            if (q > 0) { ## Simulate VARMA(p,q)
                extend.psi <- array(c(diag(k), -theta, rep(0, k * k *(n - q))), dim = c(k, k, n + 1))
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
                     temp2 <- temp2 + crossprod(t(phi[, , j]),FirstSim.Series[i - j, ])
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