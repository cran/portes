"MahdiMcLeod" <-
function(obj,lags=seq(5,30,5),fitdf=0,sqrd.res=FALSE){

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
     k <- NCOL(res)
     n <- NROW(res)

    Det <- numeric(length(lags))
    mat <- ToeplitzBlock(res, Maxlag=max(lags))
    for (i in 1:length(lags))
    Det[i] <- (-3*n/(2*lags[i]+1))*log(det(mat[(1:((lags[i] +1 ) * k)), (1:((lags[i] + 1) * k))]))
    df <- k^2*(1.5*lags*(lags+1)/(2*lags+1)-fitdf)
    NegativeDF <- which(df<0)
    df[NegativeDF] <- 0
    PVAL <- 1 - stats::pchisq(Det,df)
    PVAL[NegativeDF] <- NA
    summary <- matrix(c(lags,Det,df,PVAL),ncol=4)
    dimnames(summary) <- list(rep("", length(Det)),c("lags","statistic","df","p-value"))
  return(summary)
}
