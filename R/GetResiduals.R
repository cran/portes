"GetResiduals" <-
  function(obj)
{
    if (class(obj) != "ar" && class(obj) != "arima0" && class(obj) != "Arima" && class(obj) != "varest" && class(obj) != "FitAR" && class(obj) != "FitFGN") 
        stop("must be class ar, arima0, Arima, varest, FitAR, or FitFGN object")
    if (class(obj)=="ar"){
            order <- obj$order
        res <- ts(as.matrix(obj$resid)[-(1:order),])
    }
      else if (all(class(obj) == "arima0") || all(class(obj) == "Arima")) {
	  pdq <- obj$arma
	  p <- pdq[1]
	  q <- pdq[2]
	  d <- pdq[6]
        order <- p+q
         res <- ts(obj$residuals) 
    }
    else if (class(obj)=="varest"){
     order <- obj$p
     res <- resid(obj)
    }
    else if (class(obj)=="FitAR"){ 
     order <- length(obj$phiHat)
     res <- ts(obj$res) 
    }
    else if (class(obj) == "FitFGN") {
        order <- 0
        res <- ts(obj$res)
    }
  return(list(order = order, res = res))
}

