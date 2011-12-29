\name{destFromBearing}
\alias{destFromBearing}
\title{
Find destination from point of origin given bearing and distance
}
\description{
Find destination from point of origin given bearing and distance
}
\usage{
destFromBearing(lon, lat, bearing, distance)
}
\arguments{
  \item{lon}{
Longitude of origin
}
  \item{lat}{
Latitude of origin
}
  \item{bearing}{
Bearing to destination
}
  \item{distance}{
Distance to cover to destination in km
}
}
\value{
Returnes the destination point(s) as matrix}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{
\code{\link{addWidth}}
}
\examples{
res <- destFromBearing(rep(2.5,10),rep(51.5,10),
                       runif(10,0,360),runif(10,0,0.1))
plot(res[,1],res[,2])
points(2.5,51.5,cex=1.1,pch=19)
}
