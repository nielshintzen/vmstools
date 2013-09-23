\name{overlapPolygons}
\alias{overlapPolygons}
\title{
Calculate the surface area that polygons have in common
}
\description{
Calculate the surface area that 2 or more polygons have in common
}
\usage{
overlapPolygons(pol1 = NULL, pol2 = NULL, projection = "LL", zone = NULL)
}
\arguments{
  \item{pol1}{
Polygon number 1. Can be of class 'data.frame' with first column Longitude and
second column Latitude. Can be of class 'PolySet' and can be of class 'SpatialPolygons'.
}
  \item{pol2}{
Polygon number 2. Can be of class 'data.frame' with first column Longitude and
second column Latitude. Can be of class 'PolySet' and can be of class 'SpatialPolygons'.
}
  \item{projection}{
Optional projection attribute to add (often "LL" for Longitude-Latitude).
}
  \item{zone}{
Optional zone attribute to add.
}
}
\value{
Returns a data.frame with overlap in km2 (or non-projected units) with PID
referring to the combination of the two polygon sets.
}
\author{
Katell Hamon, Niels Hintzen
}
\seealso{
\code{\link{lonLat2SpatialPolygons}}, \code{\link{as.PolySet}}, \code{\link{surface}}
}
\examples{


#- Test with data.frame polygons
pol1 <- data.frame(cbind(c(2,3,3,2),c(54,54,55,55)))
pol2 <- data.frame(cbind(c(2,3,3,2),c(55,55,54.5,54.5)))

overlapPolygons(pol1,pol2)

#- Test with SpatialPolygons
pol1 <- lonLat2SpatialPolygons(SI_LONG=c(2,3,3,2),SI_LATI=c(54,54,55,55))
pol2 <- lonLat2SpatialPolygons(SI_LONG=c(2,3,3,2),SI_LATI=c(54.5,54.5,55,55))

overlapPolygons(pol1,pol2)
surface(pol1)@polygons[[1]]@Polygons[[1]]@area

#- Test with PolySet polygons
pol1 <- as.PolySet(data.frame(PID=rep(1,4),POS=1:4,X=c(2,3,3,2),Y=c(54,54,55,55)))
pol2 <- as.PolySet(data.frame(PID=rep(1,4),POS=1:4,X=c(2,3,3,2),Y=c(54.5,54.5,55,55)))

overlapPolygons(pol1,pol2)

#- Test with multiple polygons
data(tacsat)
pols1 <- cbind(s1=rep(2,length(seq(49,63,2))),s2=c(seq(49,63,2)))
pols2 <- cbind(s1=tacsat$SI_LONG[seq(2,nrow(tacsat),length.out=5)],s2=tacsat$SI_LATI[seq(2,nrow(tacsat),length.out=5)])
resx   <- 1; resy <- 0.5
sPols1                       <- lonLat2SpatialPolygons(lst=lapply(as.list(1:nrow(pols1)),
                                  function(x){data.frame(SI_LONG=c(pols1[x,"s1"]-resx/2,rep(pols1[x,"s1"]+resx/2,2),pols1[x,"s1"]-resx/2),
                                                         SI_LATI=c(rep(pols1[x,"s2"]-resy/2,2),rep(pols1[x,"s2"]+resy/2,2)))}))
sPols2                       <- lonLat2SpatialPolygons(lst=lapply(as.list(1:nrow(pols2)),
                                  function(x){data.frame(SI_LONG=c(pols2[x,"s1"]-resx/2,rep(pols2[x,"s1"]+resx/2,2),pols2[x,"s1"]-resx/2),
                                                         SI_LATI=c(rep(pols2[x,"s2"]-resy/2,2),rep(pols2[x,"s2"]+resy/2,2)))}))

overlapPolygons(sPols1,sPols2)

}