"InvertQ" <- 
function(coef){
 stopifnot(is.array(coef)||dim(coef)[1] == dim(coef)[2])
 k <- dim(coef)[1]
 order <- dim(coef)[3]
  if (order==1){
   ans <- eigen(coef[,,1], symmetric=FALSE, only.values =TRUE)$value 
    MaxEigenvalue <- max(Mod(ans))
     if (MaxEigenvalue >= 1) 
      return( warning("check stationary/invertibility condition !"))
  }
 else{
    out <- matrix(numeric((k*order)^2),k*order,k*order)
      out[1:k,] <- coef
      mat <- diag((order-1)*k)
      out[((k+1):(k*order)),1:((order-1)*k)] <- mat
      ans <- eigen(out, symmetric=FALSE, only.values =TRUE)$value 
      MaxEigenvalue <- max(Mod(ans))
      if (MaxEigenvalue >= 1) 
      return( warning("check stationary/invertibility condition !"))
  }
}