"ImpulseVMA" <-
function(ar=NULL,ma=NULL,trunc.lag=NULL){
    if (!is.null(ar) && !inherits(ar,"array") && !inherits(ar,"numeric"))
	stop("ar must be enterd as NULL or array with dimension (k*k*p) or numeric")
    if (!is.null(ma) && !inherits(ma,"array") && !inherits(ma,"numeric"))
	stop("ma must be enterd as NULL or array with dimension (k*k*q) or numeric")
    if (all(ar == 0))
      ar <- NULL
     if (all(ma == 0))
      ma <- NULL
    if (is.null(ar)  && is.null(ma)) 
       return(NULL)
    if (inherits(ar,"numeric"))
      ar <- array(ar,dim=c(1,1,length(ar)))
   if (inherits(ma,"numeric"))
      ma <- array(ma,dim=c(1,1,length(ma)))
    p <- ifelse(is.null(ar),0,dim(ar)[3])
    q <- ifelse(is.null(ma),0,dim(ma)[3])
    if (is.null(trunc.lag)) trunc.lag <- p + q
    if (trunc.lag < p + q)
     stop("'truncation lag' must be as long as 'P + q'")
    k <- ifelse(p > 0 ,NROW(ar[,,1]),NROW(ma[,,1]))
      if (p==0) {
        InvertQ(ma)
        psi <- array(c(diag(k),ma,rep(0,k*k*trunc.lag)),dim=c(k,k,q+trunc.lag+1))[,,1:(trunc.lag+1)]
        return(array(psi,dim=c(k,k,trunc.lag)))
      }
   if (p>0 && q==0){
      InvertQ(ar)
       psi <- array(c(diag(k),numeric(k*k*trunc.lag)), dim=c(k,k,trunc.lag+1))
       for(j in 2:(trunc.lag+1)){
         psij <- matrix(rep(0, k),k,k)
          for(i in 1:min(j-1,p))
            psij <- psij + crossprod(t(ar[,,i]), psi[,,j-i])
          psi[,,j] <- psij
       }
     return(psi)
   }
   else {
      InvertQ(ar)
      InvertQ(ma)
        psi <- array(c(diag(k),numeric(k*k*trunc.lag)), dim=c(k,k,trunc.lag+1))
        psi[,,2:(q+1)] <- ma
      for(j in 2:(trunc.lag+1)){
       psij <- matrix(rep(0, k),k,k)
        for(i in 1:min(j-1,p))
         psij <- psij + crossprod(t(ar[,,i]), psi[,,j-i])
        psi[,,j] <- psij+psi[,,j]
      }
     return(psi)
    }
}