\name{maxRangeCI}
\alias{maxRangeCI}
\title{Computation of outer range of vessel trajectory between two succeeding points}
\description{
A vessel can, based on an assumed speed and time interval, only travel a certain maximum
range between two succeeding VMS / GPS points. This outer range can be describes as an
ellipse surrounding these two points. The outer region returned is a matrix of 360 x-y coordinates.
As well, the maximum distance able to travel depending on the speed and time interval is returned as
is a warning if the measured distance between the two succeeding points is larger
than the possible distance able to travel in the amount of time and speed given, indicating
that the vessel must have travelled faster than the given speeds.
}
\usage{
maxRangeCI(x, y, Time, speed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Longitudes of the GPS positions}
  \item{y}{Latitudes of the GPS positions}
  \item{Time}{Time in minutes between the succeeding datapoints}
  \item{speed}{Speeds at the GPS positions in nots}
}
\details{
A list is returned with in the first element the 360 x-y coordinates of the outer region.
In the second element, the maximum distance able to travel given the speed and time interval, is returned.
The third element holds a warning if the maximum distance computed is exceeded. If warn equals 0, the maximum distance
is not exceeded. If warn equals 1, this distance is exceeded.
}
\references{Pfoser and Jensen 1999 Capturing the Uncertainty of Moving-Object Representations,
Mills et al. 2006 Estimating high resolution trawl fishing effort from satellite-based vessel monitoring system data,
Hintzen et al. 2010 Improved estimation of trawling tracks using cubic Hermite spline interpolationof position registration data}
\author{Niels T. Hintzen}
\note{Used within the plotCIinterpolation() function.
}
\seealso{\code{\link{plotCIinterpolation}}, \code{\link{interpolateTacsat}},
         \code{\link{N1p0}}, \code{\link{plotInterpolation}}}
\examples{
data(tacsat)
tacsat <- formatTacsat(tacsat)
#shorten the tacsat file for example use and add date-time combination
tacsat <- tacsat[1:10,]
lon <- tacsat$SI_LONG[c(3,4)]
lat <- tacsat$SI_LATI[c(3,4)]
#Calculate maximum range to be used within Confidence Interval calculation
timeDiff <- sum(an(unlist(strsplit(tacsat$SI_TIME[4],":"))) * c(60,1,1/60))
            - sum(an(unlist(strsplit(tacsat$SI_TIME[3],":"))) * c(60,1,1/60))
res      <- maxRangeCI(lon,lat,timeDiff,tacsat$SI_SP[c(3,4)])

}

