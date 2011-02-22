\name{surface}
\alias{surface}
\title{
Calculate surface from grid cells
}
\description{
Calculate the surface in km2 of the grid that has been used
}
\usage{
surface(vmsGrid, method = "Trapezoid", includeNA = T)
}
\arguments{
  \item{vmsGrid}{defined SpatialPixelsDataFrame object}
  \item{method}{Method to be used to calculate surface, either Trapezoid or UTM}
  \item{includeNA}{Whether to include cells which do not hold any data}
}
\details{
Method UTM might take longer due to way of calculation, but is more precise than the Trapezoid function, especially when larger gridcells are used
}
\value{
The vmsGrid is returned with attached a column on cellArea
}
\references{EU lot 2 project}
\author{Fabrizio Manco, Niels T. Hintzen}
\seealso{\code{createGrid()}}
\examples{
data(tacsat)

  #Sort the Tacsat data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Get the ranges of the tacsat data
xrange  <- range(tacsat$SI_LONG,na.rm=T);
xrange  <- c(min(xrange) - min(xrange)*0.05,
             max(xrange) + max(xrange)*0.05)
yrange  <- range(tacsat$SI_LATI,na.rm=T);
yrange  <- c(min(yrange) - min(yrange)*0.05,
             max(yrange) + max(yrange)*0.05)
  #Setup a grid
grid          <- createGrid(xrange,yrange,0.1,0.05)
spatialGrid   <- SpatialGrid(grid=grid);
gridded(spatialGrid) = TRUE
sP            <- as(spatialGrid,"SpatialPixels");
sPDF          <- as(sP,"SpatialGridDataFrame")
sPDF@data     <- data.frame(rep(0,length(sPDF@grid.index)))

  #Calculate the cell surface
result        <- surface(sPDF,method="Trapezoid",includeNA=T)
print(result@data)
}
