\name{degree2Km}
\alias{degree2Km}
\title{Compute distance from degrees into kilometers}
\description{
Function transformes the distance expressed in degrees into kilometers. This based on the
GPS location of a point.
}
\usage{
degree2Km(lon, lat, degree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lon}{Longitude of the GPS position}
  \item{lat}{Latitude of the GPS positiona}
  \item{degree}{Value in degrees to turn into Km}
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\note{Computation of Km is approximation based on the Haversine formula
}
\seealso{\code{\link{distance}}, \code{\link{km2Degree}}, \code{\link{lonLatRatio}}}
\examples{
lon <- -4
lat <- 50
degree <- 1.601833

degree2Km(lon,lat,degree) #114.4897km
}

