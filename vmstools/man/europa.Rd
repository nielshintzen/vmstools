\name{europa}
\alias{europa}
\docType{data}
\title{
Polygon dataset of Europe}
\description{
Polygon dataset with the coastlines of Europe.
}
\usage{data(europa)}
\format{
  A data frame with 83878 observations on the following 5 variables.
  \describe{
    \item{\code{PID}}{a numeric vector}
    \item{\code{SID}}{a numeric vector}
    \item{\code{POS}}{a numeric vector}
    \item{\code{X}}{Longitude values}
    \item{\code{Y}}{Latitude values}
  }
}
\source{
NOAA coastline datafile}
\examples{
data(europa)
require(PBSmapping)

eurPols     <- lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),
                        function(x){data.frame(
                          SI_LONG=subset(europa,SID==x)$X,
                          SI_LATI=subset(europa,SID==x)$Y)}))

plot(eurPols,col="green",xlim=c(-4,10),ylim=c(48,62))
}
\keyword{datasets}
