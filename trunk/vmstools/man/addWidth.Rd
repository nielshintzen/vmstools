\name{addWidth}
\alias{addWidth}
\title{
Add width of gear to interpolated gear track
}
\description{
Creates a SpatialPolygon object containing interpolated tracks including the
width of a gear and can be plotted as such.
}
\usage{
addWidth(interpolation, gearWidth)
}
\arguments{
  \item{interpolation}{
interpolated dataset as obtained from the function 'interpolateTacsat'
}
  \item{gearWidth}{
Width of the gear in km
}
}

\value{
Returnes an object of class 'SpatialPolygon' from the 'sp' package. Contains all
polygons per interpolation, where each interpolation consists of several polygons.
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{
\code{\link{interpolateTacsat}}, \code{\link{destFromBearing}}
}
\examples{
require(sp)

data(tacsat)

  #Sort the VMS data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Filter the Tacsat data
tacsat     <- filterTacsat(tacsat,c(4,8),hd=NULL,remDup=TRUE)

  #Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,res=100,
                  method="cHs",params=list(fm=0.5,distscale=20,sigline=0.2,
                  st=c(2,6)),headingAdjustment=0)

res           <- addWidth(interpolation,0.024)
plot(res)
}