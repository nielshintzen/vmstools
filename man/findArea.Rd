\name{findArea}
\alias{findArea}
\title{
Find area surface smaller than threshold}
\description{
Find gridcells with maximum surface equal or smaller to 'threshold' and return their
total value and gridcells involved}
\usage{
findArea(grid, threshold = 100, diagonal = TRUE)
}
\arguments{
  \item{grid}{
A SpatialGridDataFrame
}
  \item{threshold}{
Maximum value of surfaces added up
}
  \item{diagonal}{
Allow diagonal steps in gridcell selection
}
}
\details{
Diagonal means if area may include x+1,y+1 steps v versus only x+1 or y+1 steps}
\value{
data.frame with minimum surface area and gridcells selected}
\author{
Niels T. Hintzen}
\seealso{
\code{\link{surface}}}

\examples{
xrange  <- c(0,4)
yrange  <- c(52,56)
resx    <- 0.25
resy    <- 0.125

#-create grid and assign value column
grd <- createGrid(xrange,yrange,resx,resy,type="SpatialGridDataFrame",exactBorder=TRUE)
grd@data$value <- runif(nrow(coordinates(grd)),5,10)

#- find gridcells with maximum surface equal or smaller to 'threshold' and return their
#  total value and gridcells involved (diagonal means if area may include x+1,y+1 steps v
#  versus only x+1 or y+1 steps). Note, the larger the threshold, the longer the function
#  will run!
res     <- findArea(grd,threshold=1000,diagonal=TRUE)

#- Plot the result
plot(grd,type="p",pch=19,cex=0.5)
map.axes()
selec   <- which.min(res$minval)
points(coordinates(grd)[as.numeric(unlist(strsplit(res[selec,"idxs"]," "))),],col=2,lwd=3)
}