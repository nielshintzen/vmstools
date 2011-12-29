\name{calculateCI}
\alias{calculateCI}
\title{Calculate the CI between two succeeding VMS datapoints}
\description{
The interpolated tracks can be surrounded by a sort of confidence interval representing
the outer region a vessel could have travelled between two succeeding datapoints. Within
this function the CI's are computed.

}
\usage{    
calculateCI(intLon,intLat,vmsIdx1,vmsIdx2,tacsat,
            grid,sPDF,interpolation,int,params)
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
  \item{params}{list of parameters used to perform interpolation}
}
\value{
Returns a list object with 8 elements:
1- Values of CI
2- Identifyer of grid cells that have been used in calculation
3- Maximum range vessel can have travelled given speed and time
4- grid as given in input or changed if grid needed to be enlarged to fit the calculation
5- SpatialGridDataFrame as given in input or changed if needed to be enlarged
6- SpatialGrid
7- Shortest distance of each grid cell midpoint to interpolated line
8- Shortest distance of each grid cell midpoint to VMS pings
}
\references{Hintzen et al. 2010 Improved estimation of trawling tracks using cubic Hermite spline interpolation of position registration data,
EU lot 2 project}
\author{Niels T. Hintzen}
\note{Computation can take a very long time if either grid resolution is high or
if many interpolations are used.
}
\seealso{\code{\link{interpolateTacsat}}}
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


res <- calculateCI(intLon=lon,intLat=lat,vmsIdx1=12,vmsIdx2=13,VMS.=tacsat,
                   grid=grid,sPDF=sPDF,interpolation=interpolation,int=4,
                   params=list(fm=0.45,distscale=3.1,sigline=0.42,st=c(2,6)))
image(matrix(res[[1]],nrow=67,ncol=54))
}
