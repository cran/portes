"GVStat" <-
function(obj,lags=seq(5,30,5),order=0,SquaredQ=FALSE){
     TestType<-"0"
     if (class(obj)=="ts"||class(obj)=="numeric"||class(obj)=="matrix"||(class(obj)[1]=="mts" && class(obj)[2]=="ts"))
           TestType<-"1"
     if (class(obj) == "ar" || class(obj) == "Arima" || class(obj) == "arima0" || class(obj) == "FitAR" || class(obj) == "FitFGN") 
           TestType<-"2"
     if (TestType == "0")
           stop("obj must be class ar, Arima, arima0, FitAR, FitFGN, ts, numeric, matrix, or (mts ts)")
     Maxlag <- max(lags)
     if (TestType=="1")
       res <- as.ts(obj)
     else{
             GetResid <- Get.Resid(obj)
             res <- GetResid$res
             order <- GetResid$order
       }
     if (SquaredQ) res <- res^2
     k <- NCOL(res)
     n <- NROW(res)
    Det <- numeric(length(lags))
    mat <- blockToeplitz(res, lag.max=max(lags))
    for (i in 1:length(lags)) 
     Det[i] <- (-2*n/(lags[i]+1))*log(det(mat[(1:((lags[i] +1 ) * k)), (1:((lags[i] + 1) * k))]))
    df <- k^2*(lags-order)
    df[which(df<0)] <- 0
    PVAL <- 1 - pchisq(Det,df)
    summary <- matrix(c(lags,Det,df,PVAL),ncol=4)
    dimnames(summary) <- list(rep("", length(Det)),c("Lags","Statistic","df","p-value"))
  return(summary)
}