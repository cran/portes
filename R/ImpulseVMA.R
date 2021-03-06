"ImpulseVMA" <-
function(phi=NULL,theta=NULL,trunc.lag=NULL){
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
    if (is.null(trunc.lag)) trunc.lag <- p + q
    if (trunc.lag < p + q)
     stop("'truncation lag' must be as long as 'P + q'")
    k <- ifelse(p > 0 ,NROW(phi[,,1]),NROW(theta[,,1]))
      if (p==0) {
        InvertQ(theta)
        psi <- array(c(diag(k),-theta,rep(0,k*k*trunc.lag)),dim=c(k,k,q+trunc.lag+1))[,,1:(trunc.lag+1)]
        return(array(psi,dim=c(k,k,trunc.lag)))
      }
   if (p>0 && q==0){
      InvertQ(phi)
       psi <- array(c(diag(k),numeric(k*k*trunc.lag)), dim=c(k,k,trunc.lag+1))
       for(j in 2:(trunc.lag+1)){
         psij <- matrix(rep(0, k),k,k)
          for(i in 1:min(j-1,p))
            psij <- psij + crossprod(t(phi[,,i]), psi[,,j-i])
          psi[,,j] <- psij
       }
     return(psi)
   }
   else {
      InvertQ(phi)
      InvertQ(theta)
        psi <- array(c(diag(k),numeric(k*k*trunc.lag)), dim=c(k,k,trunc.lag+1))
        psi[,,2:(q+1)] <- -theta
      for(j in 2:(trunc.lag+1)){
       psij <- matrix(rep(0, k),k,k)
        for(i in 1:min(j-1,p))
         psij <- psij + crossprod(t(phi[,,i]), psi[,,j-i])
        psi[,,j] <- psij+psi[,,j]
      }
     return(psi)
    }
}