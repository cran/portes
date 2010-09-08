"Get.Resid" <-
  function(obj)
{
    if (class(obj) != "ar" && class(obj) != "Arima" && class(obj) != "arima0" && class(obj) != "FitAR" && class(obj) != "FitFGN") 
        stop("must be class ar, Arima, arima0, FitAR, or FitFGN object")
    if (class(obj)=="ar"){
            order <- obj$order
        res <- ts(as.matrix(obj$resid)[-(1:order),])
    }
    else if (class(obj) == "Arima" || class(obj) == "arima0") {
     pq<-eval(parse(text=(as.character(obj$call))[[3]]))
      p<-pq[1]
      q<-pq[3]
      order <- p+q
      res <- ts(obj$residuals) 
    } 
    else if (class(obj)=="FitAR"){ 
     order <- length(obj$phiHat)
     res <- ts(obj$res) 
    }
    else if (class(obj) == "FitFGN") {
        order <- 1
        res <- ts(obj$res)
    }
  return(list(order = order, res = res))
}

