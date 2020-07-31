\name{lonLat2SpatialPolygons}
\alias{lonLat2SpatialPolygons}
\title{
Creates a 'SpatialPolygons' object from longitude-latitude dataset}
\description{
Takes one set of longitude-lattitude or a list of lon-lat objects and
converts it into one SpatialPolygons object.
}
\usage{
lonLat2SpatialPolygons(SI_LONG, SI_LATI, lst)
}
\arguments{
  \item{SI_LONG}{
Set with longitude values}
  \item{SI_LATI}{
Set with latitude values
}
  \item{lst}{
list of objects, each with longitude and latitude values
}
}
\details{
If the longitude and latitude set do not make a full round, the function
will do that for you.
}
\value{
Returns a 'SpatialPolygons' set from the package 'sp'}
\author{
Niels T. Hintzen}
\note{
May take a while when many polygons need to be created
}
\seealso{
\code{\link{SpatialPolygons}}
}
\examples{
data(europa)

eurPols     <- lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),
                        function(x){data.frame(
                          SI_LONG=subset(europa,SID==x)$X,
                          SI_LATI=subset(europa,SID==x)$Y)}))

area        <- lonLat2SpatialPolygons(SI_LONG=c(2,2.5,2.7,2.1),SI_LATI=c(54,54.2,55.8,55.6))
plot(eurPols,col="green",xlim=c(-4,10),ylim=c(48,62)); axis(1);axis(2);box()
plot(area,add=TRUE,col="red")
}