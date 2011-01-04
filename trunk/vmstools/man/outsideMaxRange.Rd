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
outsideMaxRange(intLon, intLat, vmsIdx1, vmsIdx2, VMS, grid, sPDF, interpolation, int, params)
}
%- maybe also 'usage' for other objects documented here.
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
\seealso{\code{createGrid(),calculateCI(),point.in.polygon()}}
\examples{
data(tacsat)

  #Sort the Tacsat data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Filter the Tacsat data
tacsat          <- filterTacsat(tacsat,c(2,6),hd=NULL,remDup=T)
tacsat$SI_DATIM <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M:%S")

  #Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,res=100,method="cHs",params=list(fm=0.5,distscale=20,sigline=0.2,st=c(2,6)),headingAdjustment=0)

  #Create the final grid where all interpolations should fit on
xrange  <- range(unlist(lapply(lapply(interpolation,function(x){return(x[-1,])}),function(x){return(x[,1])})),na.rm=T)
yrange  <- range(unlist(lapply(lapply(interpolation,function(x){return(x[-1,])}),function(x){return(x[,2])})),na.rm=T)
xrange  <- range(c(xrange,range(tacsat$SI_LONG,na.rm=T))); xrange <- c(min(xrange) - min(xrange)*0.05,
                                                                       max(xrange) + max(xrange)*0.05)
yrange  <- range(c(yrange,range(tacsat$SI_LATI,na.rm=T))); yrange <- c(min(yrange) - min(yrange)*0.05,
                                                                       max(yrange) + max(yrange)*0.05)
grid          <- createGrid(xrange,yrange,0.1,0.05)
spatialGrid   <- SpatialGrid(grid=grid);
gridded(spatialGrid) = TRUE
sP            <- as(spatialGrid,"SpatialPixels");
sPDF          <- as(sP,"SpatialPixelsDataFrame")
sPDF@data     <- data.frame(rep(0,length(sPDF@grid.index)))

lon           <- tacsat$SI_LONG[c(3,4)]
lat           <- tacsat$SI_LATI[c(3,4)]

result <- outsideMaxRange(intLon=lon,intLat=lat,vmsIdx1=3,vmsIdx2=4,tacsat,grid,sPDF,interpolation,10,params=list(fm=0.5,distscale=20,sigline=0.2,st=c(2,6)))
print(result[[2]]/result[[1]]) #0% of the total value of the CI is located outside the maximum range

}
