\name{calculateCI}
\alias{calculateCI}
\title{Calculate the CI between two succeeding VMS datapoints}
\description{
The interpolated tracks can be surrounded by a sort of confidence interval representing
the outer region a vessel could have travelled between two succeeding datapoints. Within
this function the CI's are computed.

}
\usage{    
calculateCI(intLon,intLat,vmsIdx1,vmsIdx2,tacsat,grid,sPDF,interpolation,int,params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intLon}{set of two succeeding longitude points used within the calculation}
  \item{intLat}{set of two succeeding longitude points used within the calculation}
  \item{vmsIdx1}{row number referring to the first GPS points as given in intLon/intLat}
  \item{vmsIdx2}{row number referring to the second GPS points as given in intLon/intLat}
  \item{tacsat}{tacsat dataset}
  \item{grid}{Specified grid as gridTopology (see also: 'createGrid()')}
  \item{sPDF}{Spatial dataframe to be used to allign to grid}
  \item{interpolation}{interpolation dataset}
  \item{int}{list index of interpolation to use}
  \item(params}{list of parameters used to perform interpolation}
}
\references{Hintzen et al. 2010 Improved estimation of trawling tracks using cubic Hermite spline interpolation of position registration data,
EU lot 2 project}
\author{Niels T. Hintzen}
\note{Computation can take a very long time if either grid resolution is high or
if many interpolations are used.
}
\seealso{\code{interpolateTacsat()}}
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


lon <- tacsat$SI_LONG[c(3,4)]
lat <- tacsat$SI_LATI[c(3,4)]


calculateCI(lon,lat,3,4,tacsat,grid,sPDF,interpolation,3,list(fm=0.5,distscale=20,sigline=0.2,st=c(2,6)))
}
