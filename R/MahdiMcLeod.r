"MahdiMcLeod" <-
function(obj,lags=seq(5,30,5),order=0,season=1,squared.residuals=FALSE){
     TestType <- "0"
    if (class(obj) == "ts" || class(obj) == "numeric" || class(obj) == 
        "matrix" || (class(obj)[1] == "mts" && class(obj)[2] == "ts")) 
        TestType <- "1"
    if (class(obj) == "ar" || class(obj) == "arima0" || class(obj) == 
        "Arima" || (class(obj)[1] == "ARIMA" && class(obj)[2] == "Arima") || class(obj) == "varest" || class(obj) == "lm"
        || (class(obj)[1] == "glm" && class(obj)[2] == "lm") || class(obj) == "list") 
        TestType<-"2"
    if (TestType == "0") 
        stop("obj must be class ar, arima0, Arima, (ARIMA Arima), varest, lm, (glm lm), ts, numeric, matrix, (mts ts), or list")
    Maxlag <- max(lags)
     if (TestType=="1")
       res <- as.ts(obj)
     else{ 
          GetResid <- GetResiduals(obj)
          res <- GetResid$res
          order <- GetResid$order
     }
     if (squared.residuals) 
         res <- res^2
     k <- NCOL(res)
     n <- NROW(res)
       if (Maxlag*season >= n)
         stop("Maximum value of arguments lags * season can't exceed n!")
    Det <- numeric(length(lags))
    mat <- ToeplitzBlock(res, lag.max=max(lags),season=season)
    for (i in 1:length(lags))
    Det[i] <- (-3*n/(2*lags[i]+1))*log(det(mat[(1:((lags[i] +1 ) * k)), (1:((lags[i] + 1) * k))]))
    df <- k^2*(1.5*lags*(lags+1)/(2*lags+1)-order)
    NegativeDF <- which(df<0)
    df[NegativeDF] <- 0
    PVAL <- 1 - stats::pchisq(Det,df)
    PVAL[NegativeDF] <- NA
    summary <- matrix(c(lags,Det,df,PVAL),ncol=4)
    dimnames(summary) <- list(rep("", length(Det)),c("lags","statistic","df","p-value"))
  return(summary)
}