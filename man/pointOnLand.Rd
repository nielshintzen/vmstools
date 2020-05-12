\name{pointOnLand}
\alias{pointOnLand}
\title{
Find points on land given a set of coordinates}
\description{
Find the points that are on land given a set of coordinates and polygons that
determine the part that is land
}
\usage{
pointOnLand(tacsat, lands,proj4string=NULL)
}
\arguments{
  \item{tacsat}{
Tacsat file
}
  \item{lands}{
Polygon of area that is considered to be land
}
  \item{proj4string}{Projection string, default to NULL.}
}
\details{
With many coordinates, the checking might take longer.

}
\value{
Returns a vector with values 0 and 1. 1 indicating points on land, 0 indicating points not on land.
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{
\code{\link{pointInHarbour}}
}
\examples{

data(tacsat)
data(europa)
tacsat  <- tacsat[1:1000,]
tacsat  <- sortTacsat(tacsat)

pols    <- lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),
                function(x){data.frame(SI_LONG=subset(europa,SID==x)$X,
                                       SI_LATI=subset(europa,SID==x)$Y)}))
idx     <- pointOnLand(tacsat,pols);
idx     <- which(idx == 1)

plotMap(europa,xlim=c(0,10),ylim=c(48,62))
points(tacsat$SI_LONG[idx],tacsat$SI_LATI[idx],col="red",cex=0.5,pch=19)
}