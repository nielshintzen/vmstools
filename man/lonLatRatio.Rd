\name{lonLatRatio}
\alias{lonLatRatio}
\title{Calculate the ratio between 1 degree in longitude versus 1 degree in latitude}
\description{
The distance in Km on the longitude direction changes along the latitude direction.
This function computes the ratio between 1 degree in the longitude direction depending
on the latitude of the GPS position. Returns the ratio's of two GPS locations (two
succeeding VMS datapoints). Can be used with 1 GPS position too, return NA for second value.
}
\usage{
lonLatRatio(lon, lat)
}
\arguments{
  \item{lon}{Longitude of the two GPS positions}
  \item{lat}{Latitude of the two GPS positions}
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\note{Computation is approximation based on the Haversine formula
}
\seealso{\code{\link{distance}}, \code{\link{degree2Km}}, \code{\link{km2Degree}}}
\examples{
lon <- -4
lat <- 50

lonLatRatio(lon,lat)
}

