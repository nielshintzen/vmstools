\name{CSquare}
\alias{CSquare}
\title{Calculate the CSquare notation from the GPS locations}
\description{
Compute the CSquare notation from the GPS location where you can define what resolution
you want the CSquare notation to be.
}
\usage{
CSquare(lon, lat, degrees)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lon}{Longitudes of points}
  \item{lat}{Latitudes of points}
  \item{degrees}{Resolution of CSquare notation: 10, 5, 1, 0.5, 0.1, 0.05, 0.01}
}
\references{EU Lot 2 project, based on CSquare: http://www.marine.csiro.au/csquares/ }
\author{Niels T. Hintzen}
\seealso{\code{\link{km2Degree}}, \code{\link{degree2Km}}, \code{\link{lonLatRatio}},  \code{\link{distance}}}
\examples{
lon <- -4
lat <- 50
degrees <- 0.01

CSquare(lon,lat,degrees) # "7500:104:100:100"
}
