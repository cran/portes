"BoxPierce" <-
function(obj,lags=seq(5,30,5),order=0,SquaredQ=FALSE){
     TestType<-"0"
     if (class(obj)=="ts"||class(obj)=="numeric"||class(obj)=="matrix"||(class(obj)[1]=="mts" && class(obj)[2]=="ts"))
           TestType<-"1"
     if (class(obj) == "ar" || class(obj) == "arima0" || class(obj) == "Arima" || class(obj) == "varest" || class(obj) == "FitAR" || class(obj) == "FitFGN") 
           TestType<-"2"
    if (TestType == "0")
           stop("obj must be class ar, arima0, Arima, varest, FitAR, FitFGN, ts, numeric, matrix, or (mts ts)")
    Maxlag <- max(lags)
    if (TestType=="1")
       res <- as.ts(obj)
     else{
             GetResid <- GetResiduals(obj)
             res <- GetResid$res
             order <- GetResid$order
      }
     if (SquaredQ) res <- res^2
       n <- NROW(res)
       k <- NCOL(res)
        df <- k^2*(lags-order)
        NegativeDF <- which(df<0)
        df[NegativeDF] <- 0
  	     Accmat <- acf(res, lag.max = Maxlag, plot = FALSE, type = "correlation")$acf
	     inveseR0 <- solve(Accmat[1,,])
             prodvec <- numeric(Maxlag)
        for(l in 1:Maxlag){
            tvecR <- t(as.vector(Accmat[l+1,,]))
	        prodvec[l] <- crossprod(t(tvecR),crossprod(t(kronecker(inveseR0,inveseR0)),t(tvecR)))
	}
       Q <- n*cumsum(prodvec)
       STATISTIC <- Q[lags]
       PVAL <- 1 - pchisq(STATISTIC,df)
       PVAL[NegativeDF] <- NA
       summary <- matrix(c(lags,STATISTIC,df,PVAL),ncol=4)
       dimnames(summary) <- list(rep("", length(STATISTIC)),c("Lags","Statistic","df","p-value"))
    return(summary)
}