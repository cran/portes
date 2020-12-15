"GetResiduals" <-
  function(obj)
{
    class.obj = class(obj)[1]
    if (class.obj != "ar" && class.obj != "arima0" && class.obj != "Arima" && class.obj != "varest" && 
         class.obj != "ARIMA" && class.obj != "lm" 
       && class.obj != "glm" && class.obj != "list" ) 
   stop("obj must be class ar, arima0, Arima, (ARIMA forecast_ARIMA Arima), varest, lm, (glm lm), or list")      
    if (all(class.obj=="ar")){
        order <- obj$order
        res <- ts(as.matrix(obj$resid)[-(1:order),])
    }
    else if (all(class.obj == "arima0") || all(class.obj == "Arima")|| all (class.obj == "ARIMA")) {
	  pdq <- obj$arma
	  p <- pdq[1]
	  q <- pdq[2]
	  ps <- pdq[3]
	  qs <- pdq[4]
        order <- p+q+ps+qs
         res <- ts(obj$residuals) 
    }
    else if (all(class.obj=="varest")){
     order <- obj$p
     res <- resid(obj)
    }
    else if (all(class.obj == "list")){
        order <- obj$order
        if(is.null(order))
          order <- 0
        else 
          order <- order
        res <- obj$res
    }
    if (all(class.obj=="lm") || all(class.obj == "glm")){
     order <- 0
     res <- obj$residuals
    }
  return(list(order = order, res = res))
}