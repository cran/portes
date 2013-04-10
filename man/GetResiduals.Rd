\name{GetResiduals}
\alias{GetResiduals}

\title{Extract Residuals from Fitted ARIMA, VAR, or FGN Model }

\description{
This utility function is useful to use in the portmanteau functions, 
\code{\link{BoxPierce}}, \code{\link{gvtest}}, \code{\link{Hosking}},
\code{\link{LiMcLeod}}, \code{\link{LjungBox}}, and \code{\link{portest}}.
It takes a fitted time-series object with class \code{"ar"}, \code{"arima0"}, \code{"Arima"}, 
\code{"varest"}, \code{"FitAR"}, or \code{"FitFGN"} and returns the residuals and the order of the fitted object.
}

\usage{
GetResiduals(obj)
}

\arguments{
\item{obj}{a fitted time-series model with class \code{"ar"}, \code{"arima0"}, \code{"Arima"}, 
   \code{"varest"}, \code{"FitAR"}, or \code{"FitFGN"}. }
}

\value{
List of order of fitted time series model and residuals from this model. 
The class of the residuals is \code{\link{ts}}. 
}

\author{ Esam Mahdi and A.I. McLeod.}

\seealso{
\code{\link{BoxPierce}}, \code{\link{gvtest}}, \code{\link{Hosking}},
\code{\link{LiMcLeod}}, \code{\link{LjungBox}}, \code{\link{ar}}, \code{\link{arima0}},
\code{\link{arima}}, \code{\link[forecast]{Arima}}, \code{\link{FitAR}}, \code{VAR}, \code{FitFGN}
}

\examples{
fit <- arima(Nile, c(1, 0, 1))
GetResiduals(fit)
}

\keyword{ts}

