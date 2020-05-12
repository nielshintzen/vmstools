\name{km2Degree}
\alias{km2Degree}

\title{Compute distance from kilometers into degrees}
\description{
Function transformes the distance expressed in kilometers into degrees. This based on the
GPS location of a point.
}
\usage{
km2Degree(lon, lat, km)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lon}{Longitude of the GPS position}
  \item{lat}{Latitude of the GPS positiona}
  \item{km}{Value in Km to turn into degrees}
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\note{Computation of degrees is approximation based on the Haversine formula}
\seealso{\code{\link{distance}}, \code{\link{degree2Km}}, \code{\link{lonLatRatio}}}
\examples{
lon <- -4
lat <- 50
km  <- 114.4897

km2Degree(lon,lat,km) #1.601833
}

