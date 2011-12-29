\name{bearing}
\alias{bearing}
\title{Compute bearing between two points on a sphere (approximation of the earth) at the starting point}
\description{
Compute the bearing between two GPS locations defined in longitude and latitude notation
on the earth. The earth is assumed to have a perfect spherical shape. Bearing is
returned in compass degrees.
}
\usage{
bearing(lon, lat, lonRef, latRef)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lon}{Longitude of point 2}
  \item{lat}{Latitude of point 2}
  \item{lonRef}{Longitude of point 1}
  \item{latRef}{Latitude of point 1}
}
\references{EU Lot 2 project, based on the Haversine formula, see also: Hintzen et al. 2010 Fisheries Research }
\author{Niels T. Hintzen}
\seealso{\code{\link{km2Degree}}, \code{\link{degree2Km}}, \code{\link{lonLatRatio}},\code{\link{distance}}}
\examples{
lon <- -4
lat <- 50
lonRef <- -4.2
latRef <- 51

bearing(lon,lat,lonRef,latRef) #352.8271
}
