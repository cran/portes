"GetResiduals" <-
  function(obj)
{
    class.obj = class(obj)[1]
    if (class.obj != "ar" && class.obj != "arima0" && class.obj != "Arima" && class.obj != "varest" && 
         all(class(obj)[1] != "ARIMA" && class(obj)[2] != "Arima") && class.obj != "lm" 
       && all(class(obj)[1] != "glm" && class(obj)[2] != "lm") && class.obj != "list" ) 
   stop("obj must be class ar, arima0, Arima, (ARIMA Arima), varest, lm, (glm lm), or list")      
    if (all(class.obj=="ar")){
        order <- obj$order
        res <- ts(as.matrix(obj$resid)[-(1:order),])
    }
    else if (all(class.obj == "arima0") || all(class.obj == "Arima")|| all (class(obj)[1] == "ARIMA" && class(obj)[2] == "Arima")) {
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
    if (all(class.obj=="lm") || all(class(obj)[1] == "glm" && class(obj)[2] == "lm")){
     order <- 0
     res <- obj$residuals
    }
  return(list(order = order, res = res))
}