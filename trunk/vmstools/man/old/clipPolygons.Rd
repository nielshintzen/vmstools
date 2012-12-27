\name{clipPolygons}
\alias{clipPolygons}
\title{
Function to clip Polygons together
}
\description{
This function removes from a set of polygons A the areas overlapping with another set of polygons B.
}
\usage{
clipPolygons(shapeAll, europa)
}
\arguments{
  \item{shapeAll}{These are the polygons A (as shapefiles) that need to be clipped}
  \item{europa}{These are the polygons B (as a polyset format) used for the clipping}
}
\details{
This function is based on the joinPolys function of the PBSmapping package.
It is used in the DCF indicator 7 which calculates the total surface within a given area that is not impacted by fisheries.
A shapefile is the input to define this area (a MPA, for example). This function will remove the areas from this shapefile that might overlap with lands (stored in europa.rda), otherwise the indicator will be overestimated as land is, by definition, not impacted by fisheries.
}
\value{
A clipped resulting polyset is returned.
}
\references{EU lot 2 project}
\author{Fabrizio Manco}
\seealso{\code{indicators()}}
\examples{
\dontrun{

# Load the set of polygons used for the clipping, in our case, lands
data(europa)

shapeAll <- lonLat2SpatialPolygons(SI_LONG=c(0,1,1,0),SI_LATI=c(54,54,54.5,54.5))

# Do the clipping
clipShapeFromLand<-clipPolygons (shapeAll, europa)

# Plot the resulting clipped polygons to check
plotPolys(clipShapeFromLand)
}
}