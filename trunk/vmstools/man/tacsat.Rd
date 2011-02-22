\name{tacsat}
\alias{tacsat}
\docType{data}
\title{
Vessel Monitoring by Satellite system dataset
}
\description{
A dataset consisting of VMS registrations of disguised origin for 2 consecutice years.
}
\usage{data(tacsat)}
\format{
  A data frame with 97015 observations on the following 8 variables.
  \describe{
    \item{\code{VE_COU}}{Flag nation of vessel}
    \item{\code{VE_REF}}{Vessel reference / ID}
    \item{\code{SI_LATI}}{Latitudinal position (in decimal degrees)}
    \item{\code{SI_LONG}}{Longitudinal position (in decimal degrees)}
    \item{\code{SI_DATE}}{Date}
    \item{\code{SI_TIME}}{Time (24 hours clock)}
    \item{\code{SI_SP}}{Speed (in knots per hour)}
    \item{\code{SI_HE}}{Heading (in degrees)}
  }
}
\details{
This example dataset if for example purposes only. Without prior permission of the authors it is not allowed to use this data
other than for example non-publishable purposes.
}
\source{
NA
}
\references{
EU Lot 2 project, Niels T. Hintzen (niels.hintzen@wur.nl)
}
\examples{
data(tacsat)
tacsat <- formatTacsat(tacsat)
str(tacsat)
}
\keyword{datasets}
