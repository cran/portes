\name{GetResiduals}
\alias{GetResiduals}

\title{Extract Residuals from ARIMA, VAR, or any Simulated Fitted Time Series Model}

\description{
This utility function is useful to use in the portmanteau functions, 
\code{\link{BoxPierce}}, \code{\link{MahdiMcLeod}}, \code{\link{Hosking}},
\code{\link{LiMcLeod}}, \code{\link{LjungBox}}, and \code{\link{portest}}.
\code{GetResiduals()} function takes a fitted time-series object 
with class \code{"ar"}, \code{"arima0"}, \code{"Arima"}, 
\code{("ARIMA forecast ARIMA Arima")}, \code{"lm"}, \code{("glm" "lm")}, 
\code{"varest"}, or \code{"list"}. 
 and returns the residuals and the order from the fitted object.
}

\usage{
GetResiduals(obj)
}

\arguments{
\item{obj}{a fitted time-series model with class \code{"ar"}, \code{"arima0"}, \code{"Arima"}, 
   \code{("ARIMA forecast ARIMA Arima")}, \code{"lm"}, \code{("glm" "lm")}, 
   \code{"varest"}, or \code{"list"}. }
}

\value{
List of order of fitted time series model and residuals from this model. 
}

\author{ Esam Mahdi and A.I. McLeod.}

\seealso{
\code{\link{ar}}, \code{\link{ar.ols}}, \code{\link{ar.burg}}, 
\code{\link{ar.yw}}, \code{\link{ar.mle}}, \code{\link{arima0}}, \code{\link{arima}}, 
\code{\link[forecast]{Arima}}, \code{\link[forecast]{auto.arima}},
\code{\link[stats]{lm}}, \code{\link[stats]{glm}}, \code{\link[vars]{VAR}}, 
\code{\link{BoxPierce}}, \code{\link{LjungBox}}, \code{\link{MahdiMcLeod}}, \code{\link{Hosking}},
\code{\link{LiMcLeod}}.
}

\examples{
fit <- arima(Nile, c(1, 0, 1))
GetResiduals(fit)
}

\keyword{ts}

