\name{Get.Resid}
\alias{Get.Resid}

\title{Extract Residuals from Fitted ARIMA, VAR, or FGN Model }

\description{
This utility function is useful to use in the portmanteau functions, 
\code{\link{BoxPierce}}, \code{\link{GVStat}}, \code{\link{Hosking}},
\code{\link{LiMcLeod}}, \code{\link{LjungBox}}, and \code{\link{portes}}.
It takes a fitted time-series object with class \code{"ar"}, \code{"Arima"}, \code{"arima0"},
\code{"FitAR"}, or \code{"FitFGN"} and returns the residuals and the order of the fitted object.
}

\usage{
Get.Resid(obj)
}

\arguments{
\item{obj}{a fitted time-series model with class \code{"ar"}, \code{"Arima"}, \code{"arima0"},
   \code{"FitAR"}, or \code{"FitFGN"}. }
}

\value{
List of order of fitted time series model and residuals from this model. 
The class of the residuals is \code{\link{ts}}. 
}

\author{ Esam Mahdi & A.I. McLeod.}

\seealso{
\code{\link{BoxPierce}}, \code{\link{GVStat}}, \code{\link{Hosking}},
\code{\link{LiMcLeod}}, \code{\link{LjungBox}}, \code{\link{portes}}, \code{\link{ar}},  \code{\link{arima}},
 \code{\link{arima0}}, \code{FitAR}, \code{FitFGN}
}

\examples{
fit <- arima(Nile, c(2, 0, 0))
Get.Resid(fit)
}

\keyword{ts}

