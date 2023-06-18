"ToeplitzBlock" <- 
function(res,Maxlag)
{
     res <- as.matrix(res)
     k <- NCOL(res)
     n <- NROW(res)
     m <- Maxlag+1
     out <- matrix(numeric(k*m*k*(Maxlag+1)),nrow=k*m,ncol=k*m)
     Accmat <- stats::acf(res, lag.max = Maxlag, plot = FALSE, type = "covariance")$acf
      inveseC0 <- solve(Accmat[1,,])
      L <- t(chol(inveseC0))
       for (i in 0:Maxlag)
	  for (j in i:Maxlag){
	  out[(j*k+1):(k*(j+1)),(i*k+1):(k*(i+1))] <- crossprod(t(crossprod(L,t(Accmat[(j-i)+1,,]))),L)
          out[(i*k+1):(k*(i+1)),(j*k+1):(k*(j+1))] <- crossprod(t(crossprod(L,Accmat[(j-i)+1,,])),L)
        }
   return(out)
}