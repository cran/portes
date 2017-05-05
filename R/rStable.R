"rStable" <- 
function (n, Alpha, Beta, Scale = NULL, Location = NULL) 
{
    if (any(Alpha > 2)) 
        stop("Error: Alpha is greater than 2")
    if (any(Alpha <= 0)) 
        stop("Error: Alpha is less than or equal to 0")
    if (any(Beta < -1)) 
        stop("Error: Beta is less than -1")
    if (any(Beta > 1)) 
        stop("Error: Beta is greater than 1")
    k <- length(Alpha)
    if (length(Beta) !=k) 
        stop("Error: Alpha and Beta should have the same dimension")
    theta <- matrix( pi * (runif(k*n) - 1/2),ncol=k)
    z <- matrix(-log(runif(k*n)),ncol=k)
      result <- ans <- matrix(numeric(k*n),ncol=k)
        c <- (1 + (Beta * tan(pi * Alpha/2))^2)^(1/(2 * Alpha))
        theta0 = (1/Alpha) * atan(Beta * tan(pi * Alpha/2))
     for (i in 1:k){
       if (Alpha[i] == 1 & Beta[i] == 0) 
        result[,i] <- matrix(rcauchy(k*n),ncol=k)
       else
        result[,i] <- (c[i] * sin(Alpha[i] * (theta[,i] + theta0[i]))/(cos(theta[,i]))^(1/Alpha[i])) * (cos(theta[,i] - Alpha[i] * (theta[,i] + theta0[i]))/z[,i])^((1 - Alpha[i])/Alpha[i])
       result[,i] <- result[,i] - Beta[i] * tan(Alpha[i] * pi/2)
       if (is.null(Scale)) Scale <- rep(1,k)
       if (is.null(Location)) Location <- rep(0,k)
      ans[,i] <- result[,i] * Scale[i] + Location[i]
     }
  return(ans)
}
