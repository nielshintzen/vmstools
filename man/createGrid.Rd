\name{createGrid}
\alias{createGrid}
\title{create a Grid (GridTopology: sp) depending on x and y limits and gridcell size}
\description{
Function creates a grid with rounded (no decimal) begin and end longitude and latitude
based on the ranges given. It creates this grid with gridcell size specified.
}
\usage{
createGrid(xrange, yrange, resx, resy,type="GridTopology",exactBorder=FALSE)
}
\arguments{
  \item{xrange}{specify minimum and maximum longitude / x values which should at least be in the grid}
  \item{yrange}{specify minimum and maximum latitude / y values which should at least be in the grid}
  \item{resx}{gridcell size in degrees in the longitude / x direction}
  \item{resy}{gridcell size in degrees in the latitude / y direction}
  \item{type}{specify output class: 'GridTopology','SpatialGrid','SpatialPixels',
              'SpatialPixelsDataFrame','SpatialGridDataFrame'}
  \item{exactBorder}{Logical: if the specified resx and resy  (lower-left values)
               should be taken as exact}
}
\details{
The grid is created depending on the 'sp' package. The grid defined has default
class 'GridTopology'.
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{vmsGridCreate}}}
\examples{
x <- seq(-4.999,5,length.out=10)
y <- seq(50.002,55,length.out=10)
xrange <- range(x,na.rm=TRUE)
yrange <- range(y,na.rm=TRUE)

#returns a grid with 1001 cells in the x-direction and 101 in the y-direction
Grid <- createGrid(xrange,yrange,0.01,0.05,type="SpatialGrid")
bbox(Grid)
Grid <- createGrid(xrange,yrange,0.01,0.05,type="SpatialGrid",exactBorder=TRUE)
bbox(Grid)
}

