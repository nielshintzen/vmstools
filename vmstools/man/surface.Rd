\name{surface}
\alias{surface}
\title{
Calculate surface from grid cells or polygons
}
\description{
Calculate the surface in km2 of the grid / polygons that has been used
}
\usage{
surface(obj, method = "Trapezoid", includeNA = TRUE, zone = NULL)
}
\arguments{
  \item{obj}{defined SpatialGridDataFrame or SpatialPolygons object (see 'sp'
             package)}
  \item{method}{Method to be used to calculate surface, either Trapezoid or UTM}
  \item{includeNA}{Whether to include cells which do not hold any data}
  \item{zone}{Include UTM zone notation (or detected automatically when NULL}
}
\details{
Method UTM might take longer due to way of calculation, but is more precise than
the Trapezoid function, especially when larger gridcells are used.
}
\value{
If obj is a SpatialGridDataFrame an additional column named 'cellArea' is returned
which holds the km2 area of the grid cell.
If obj is a SpatialPolygons, each polygon area slot is filled with the area in
km2.
}
\author{Niels T. Hintzen}
\seealso{\code{\link{createGrid}}}
\examples{
data(tacsat)

  #Sort the Tacsat data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Get the ranges of the tacsat data
xrange  <- range(tacsat$SI_LONG,na.rm=TRUE);
xrange  <- c(min(xrange) - min(xrange)*0.05,
             max(xrange) + max(xrange)*0.05)
yrange  <- range(tacsat$SI_LATI,na.rm=TRUE);
yrange  <- c(min(yrange) - min(yrange)*0.05,
             max(yrange) + max(yrange)*0.05)
  #Setup a grid
sPDF    <- createGrid(xrange,yrange,resx=0.1,resy=0.05,type="SpatialGridDataFrame")

  #Setup a polygon
sP      <- lonLat2SpatialPolygons(lst=list(data.frame(
              SI_LONG=c(4,4.5,4.7,4),
              SI_LATI=c(54,54,55.5,55.7))))

  #Calculate the cell surface
result        <- surface(sPDF,method="Trapezoid",includeNA=TRUE)
print(head(result@data))
result        <- surface(sP,zone=31)
print(result@polygons[[1]]@Polygons[[1]]@area)
}
