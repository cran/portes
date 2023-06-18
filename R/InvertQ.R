"InvertQ" <- 
function(coef){
stopifnot((inherits(coef,"numeric")||inherits(coef,"matrix")||inherits(coef,"array")&&(dim(coef)[1]==dim(coef)[2])))
     if (inherits(coef,"numeric"))
       coef <- array(coef,dim=c(1,1,length(coef)))
     if (inherits(coef,"matrix"))
       coef <- array(coef,dim=c(NROW(coef),NROW(coef),1))
     k <- dim(coef)[1]
     fitdf <- dim(coef)[3]
     if (fitdf==1)
       ans <- eigen(coef[,,1], symmetric=FALSE, only.values =TRUE)$value 
     else{
        blockMat <- matrix(numeric((k*fitdf)^2),k*fitdf,k*fitdf)
        blockMat[1:k,] <- coef
        Imat <- diag((fitdf-1)*k)
        blockMat[((k+1):(k*fitdf)),1:((fitdf-1)*k)] <- Imat
        ans <- eigen(blockMat, symmetric=FALSE, only.values =TRUE)$value 
     }
   MaxEigenvalue <- max(Mod(ans))
   if (MaxEigenvalue >= 1) 
     return( warning("check stationary/invertibility condition !"))
}