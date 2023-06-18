"LjungBox" <-
function(obj,lags=seq(5,30,5),fitdf=0,sqrd.res=FALSE){
  Maxlag <- max(lags)

   class.obj = class(obj)[1]
     TestType <- "0"
    if (class.obj == "ts" || class.obj == "numeric" || class.obj == "matrix" || class.obj == "mts")
        TestType <- "1"

     if (class.obj == "ar" || class.obj == "arima0" || class.obj ==
         "Arima" || class.obj == "ARIMA" || class.obj == "forecast_ARIMA" || class.obj == "varest" || class.obj == "list")
       TestType<-"2"

     if (TestType == "0")
       stop("obj must be class ar, arima0, Arima, (forecast_ARIMA ARIMA Arima), varest, ts, numeric, matrix, (mts ts), or list")

    if (TestType=="1")
       res <- as.ts(obj)
     else{
          GetResid <- GetResiduals(obj)
          res <- GetResid$res
          fitdf <- GetResid$fitdf
     }

     if (sqrd.res)
         res <- res^2
       n <- NROW(res)
       k <- NCOL(res)

        df <- k^2*(lags-fitdf)
        NegativeDF <- which(df<0)
        df[NegativeDF] <- 0
	     Accmat <- stats::acf(res, lag.max = Maxlag, plot = FALSE, type = "correlation")$acf
	     inveseR0 <- solve(Accmat[1,,])
             prodvec <- numeric(Maxlag)
	  for(l in 1:Maxlag){
            tvecR <- t(as.vector(Accmat[l+1,,]))
	        prodvec[l] <- 1/(n-l)*crossprod(t(tvecR),crossprod(t(kronecker(inveseR0,inveseR0)),t(tvecR)))
	  }
       Q <- n*(n+2)*cumsum(prodvec)
       STATISTIC <- Q[lags]
       PVAL <- 1 - stats::pchisq(STATISTIC,df)
       PVAL[NegativeDF] <- NA
       summary <- matrix(c(lags,STATISTIC,df,PVAL),ncol=4)
       dimnames(summary) <- list(rep("", length(STATISTIC)),c("lags","statistic","df","p-value"))
    return(summary)
}
