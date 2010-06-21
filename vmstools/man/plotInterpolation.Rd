\name{plotInterpolation}
\alias{plotInterpolation}
\title{Plot the interpolated VMS dataset as lines}
\description{
The interpolated tracks are plotted directly onto a grid. If a certain grid cell is
trawled / covered / travelled more than once, values are summed.
}
\usage{
plotInterpolation(interpolation, VMS, grid, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interpolation}{List of interpolations as outputted by 'interpolateVMS()'}
  \item{VMS}{VMS dataset}
  \item{grid}{Specified grid as gridTopology (see also: 'createGrid()')}
  \item{plot}{True or False, indicating if figure should be created or if SpatialGridDataFrame should be returned}
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{plotCIinterpolation(),interpolateVMS()}}
\examples{
data(VMS)

  #Sort the VMS data
VMS     <- sortVMS(VMS)
VMS     <- VMS[1:1000,]

  #Filter the VMS data
VMS     <- filterVMS(VMS,c(2,6),hd=NULL,remDup=T)

  #Interpolate the VMS data
interpolation <- interpolateVMS(VMS,interval=120,margin=10,res=100,method="cHs",params=list(fm=0.5,distscale=20,sigline=0.2))

  #Create the final grid where all interpolations should fit on
xrange  <- range(unlist(lapply(lapply(interpolation,function(x){return(x[-1,])}),function(x){return(x[,1])})),na.rm=T)
yrange  <- range(unlist(lapply(lapply(interpolation,function(x){return(x[-1,])}),function(x){return(x[,2])})),na.rm=T)
xrange  <- range(c(xrange,range(VMS$declon,na.rm=T))); xrange <- c(min(xrange) - min(xrange)*0.05,
                                                                   max(xrange) + max(xrange)*0.05)
yrange  <- range(c(yrange,range(VMS$declat,na.rm=T))); yrange <- c(min(yrange) - min(yrange)*0.05,
                                                                   max(yrange) + max(yrange)*0.05)
grid    <- createGrid(xrange,yrange,0.1,0.05)

  #Plot the interpolated data without a CI
plotInterpolation(interpolation,VMS,grid,plot=T)
}
