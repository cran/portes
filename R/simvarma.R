"simvarma" <-
function (phi = NULL, theta = NULL, sigma, intercept = NA, n, StableParameters=NA,Trunc.Series = NA) 
{
    if (!is.null(phi) && class(phi)!="array" && class(phi)!="numeric")
	stop("Phi must be enterd as NULL or array with dimension (k*k*p) or numeric")
    if (!is.null(theta) && class(theta)!="array" && class(theta)!="numeric")
	stop("Theta must be enterd as NULL or array with dimension (k*k*q) or numeric")
    if (all(phi == 0))
      phi <- NULL
     if (all(theta == 0))
      theta <- NULL
    if (is.null(phi)  && is.null(theta)) 
       return(NULL)
    if (class(phi) == "numeric")
      phi <- array(phi,dim=c(1,1,length(phi)))
    if (class(theta)=="numeric")
      theta <- array(theta,dim=c(1,1,length(theta)))
    p <- ifelse(is.null(phi),0,dim(phi)[3])
    q <- ifelse(is.null(theta),0,dim(theta)[3])
    sigma <- as.matrix(sigma)
    k <- NCOL(sigma)
   if (p > 0 && ((dim(phi)[1] != dim(phi)[2])||dim(phi)[2] != k))
        stop("Wrong dimensions of phi or/and sigma")
   if (q > 0 && ((dim(theta)[1] != dim(theta)[2])||dim(theta)[2] != k))
        stop("Wrong dimensions of theta or/and sigma")
    if (all(is.na(intercept))) 
        intercept <- rep(0, k)
    if (length(intercept) != k) 
        stop("intercept must be entered as a vector with length equal to number of sigma rows")
    StableQ <- all(!is.na(StableParameters))
       if (StableQ) {
         StableParameters <- matrix(StableParameters,nrow=k)
           stopifnot(NCOL(StableParameters) == 4)
             ALPHA<-StableParameters[,1]
             BETA<-StableParameters[,2]
             GAMMA<-StableParameters[,3]
             DELTA<-StableParameters[,4]
       }
    if (p == 0) {
        if (StableQ)
         epsilon <- rstable(n + q, ALPHA, BETA, GAMMA, DELTA)
        else
         epsilon <- t(crossprod(chol(sigma),matrix(rnorm(k*(n + q)),ncol=n + q)))
        if (q == 0)  ## Simulate white noise
            return(ts(scale(epsilon, center=-intercept,scale = FALSE)))
         InvertQ(theta)
         ## Simulate VMA(q)
            psi <- array(c(diag(k), -theta), dim = c(k, k, q + 1))
            Sim.Series <- simvma(psi = psi, a = epsilon)
          return(ts(scale(Sim.Series, center=-intercept,scale = FALSE)))        
    }
    if (is.na(Trunc.Series)) 
        Trunc.Series <- min(100,ceiling(n/3))
    Sim.Series <- matrix(numeric(0), nrow = n, ncol = k)
    r <- max(p, q)
    psi <- ImpulseVMA(phi = phi, theta = theta, Trunc.Series = Trunc.Series)
    if (StableQ)
        epsilon <- rstable(Trunc.Series + r, ALPHA, BETA, GAMMA, DELTA)
    else
        epsilon <- t(crossprod(chol(sigma),matrix(rnorm(k*(Trunc.Series + r)),ncol=Trunc.Series + r)))
    Sim.Series[1:r, ] <- simvma(psi = psi, a = epsilon)
    a <- matrix(epsilon[1:r, ], nrow = r, ncol = k)
    if (StableQ)
       epsilon <- rbind(a,rstable(n, ALPHA, BETA, GAMMA, DELTA))
    else
        epsilon <- rbind(a,t(crossprod(chol(sigma),matrix(rnorm(k*n),ncol=n))))
     if (q > 0) { ## Simulate VARMA(p,q)
        extend.psi <- array(c(diag(k), -theta, rep(0, k * k * 
            (n - q))), dim = c(k, k, n + 1))
        u <- matrix(numeric(0), nrow = n, ncol = k)
        for (i in (q + 1):(n + q)) {
            out <- 0
            for (j in 0:q) {
                out = out + crossprod(t(extend.psi[, , j + 1]), 
                  epsilon[i - j, ])
            }
            u[i - q, ] <- out
        }
      }
    else u <- epsilon  ## Simulate VAR(p)
    for (i in (r + 1):n) {
        temp2 <- 0
        for (j in 1:p) temp2 <- temp2 + crossprod(t(phi[, , j]), 
            Sim.Series[i - j, ])
        Sim.Series[i, ] <- temp2 + u[i, ]
    }
    return(ts(scale(Sim.Series, center=-intercept,scale = FALSE))) 
}
