"GetResiduals" <-
  function(obj)
{
    class.obj = class(obj)[1]

    if (class.obj != "ar" && class.obj != "arima0" && class.obj != "Arima" && class.obj != "varest" &&
        class.obj != "ARIMA" && class.obj != "forecast_ARIMA" && class.obj != "list" )
   stop("obj must be class ar, arima0, Arima, (ARIMA forecast_ARIMA Arima), varest, or list")
    if (all(class.obj=="ar")){
        fitdf <- obj$order
        res <- ts(na.omit(obj$resid))
    }
    else if (all(class.obj == "arima0") || all(class.obj == "Arima")|| all (class.obj == "ARIMA") || all (class.obj == "forecast_ARIMA")) {
	  pdq <- obj$arma
	  p <- pdq[1]
	  q <- pdq[2]
      fitdf <- p+q
      res <- ts(obj$residuals)
    }
    else if (all(class.obj=="varest")){
     fitdf <- obj$p
     res <- resid(obj)
    }
    else if (all(class.obj == "list")){
        fitdf <- obj$fitdf
        if(is.null(fitdf))
          fitdf <- 0
        else
          fitdf <- fitdf
        res <- obj$res
    }
  return(list(fitdf = fitdf,res = res))
}
