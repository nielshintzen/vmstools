\name{outsideMaxRange}
\alias{outsideMaxRange}
\title{compute fraction of Confidence Interval that is located outside a maximum range}
\description{
The calculation of the Confidence Interval surrounding an interpolation depends on two parameters:
sigline & distscale. These use of these parameters could result in extremely wide or extremely small CI's.
To check which proportion is located inside and outside the maximum range (as defined by an ellipse), this
function calculates this proportion, as well as the maximum value representing the starting and end point
CI values.
}
\usage{
outsideMaxRange(intLon, intLat, vmsIdx1, vmsIdx2, VMS, grid,
                sPDF, interpolation, int, params)
}
\arguments{
  \item{intLon}{Longitudes of single interpolation start and endpoint}
  \item{intLat}{Latitudes of single interpolation start and endpoint}
  \item{vmsIdx1}{Start point of interpolation expressed in row number of VMS dataset}
  \item{vmsIdx2}{End point of interpolation expressed in row number of VMS dataset}
  \item{VMS}{VMS dataset}
  \item{grid}{Grid object inhereted from the 'sp' library: gridTopology}
  \item{sPDF}{Spatial Pixels Data Frame object, used to store the CI data in}
  \item{interpolation}{List of interpolations}
  \item{int}{Identifier of interpolation}
  \item{params}{parameters list(fm,distscale,sigline}
}
\details{
The spatial pixels data frame object is created from the grid object as is explained in the spatial library 'sp'.

The interpolation is the list of interpolations returned by the 'interpolateVMS' function. 
}
\references{Hintzen et al. 2010 Improved estimation of trawling tracks using cubic Hermite spline interpolation of position registration data,
EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{createGrid}}, \code{\link{calculateCI}}, \code{\link{point.in.polygon}}}
\examples{
data(tacsat)

  #Sort the Tacsat data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Filter the Tacsat data
tacsat          <- filterTacsat(tacsat,c(2,6),hd=NULL,remDup=T)

  #Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
                    res=100,method="cHs",params=list(fm=0.5,distscale=20,
                    sigline=0.2,st=c(2,6)),headingAdjustment=0)

  #Create the final grid where all interpolations should fit on
xrange        <- c(2,3); yrange <- c(51,52)
grid          <- createGrid(xrange,yrange,resx=0.01,resy=0.005)
sPDF          <- createGrid(xrange,yrange,resx=0.01,resy=0.005,type="SpatialPixelsDataFrame")
sPDF@data     <- data.frame(rep(0,length(sPDF@grid.index)))


lon <- tacsat$SI_LONG[c(12,13)]
lat <- tacsat$SI_LATI[c(12,13)]

result <- outsideMaxRange(intLon=lon,intLat=lat,vmsIdx1=12,vmsIdx2=13,
            tacsat,grid,sPDF,interpolation,4,params=
            list(fm=0.45,distscale=3.1,sigline=0.42,st=c(2,6)))

#3.7\% of the total value of the CI is located outside the maximum range
print(result[[2]]/result[[1]])

}
