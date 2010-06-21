\name{plotCIinterpolation}
\alias{plotCIinterpolation}
\title{Plot the interpolated VMS dataset including the CI's}
\description{
The interpolated tracks can be surrounded by a sort of confidence interval representing
the outer region a vessel could have travelled between two succeeding datapoints. Within
this function the CI's are computed and added together to one figure representing the average
number of times the area has possibly been trawled / covered / travelled.

}
\usage{
plotCIinterpolation(interpolation, VMS = NULL, grid, plot = F, params = list(0.2, 20, 0.2), plotPoint = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interpolation}{List of interpolations as outputted by 'interpolateVMS()'}
  \item{VMS}{VMS dataset}
  \item{grid}{Specified grid as gridTopology (see also: 'createGrid()')}
  \item{plot}{True or False, indicating if figure should be created or if SpatialGridDataFrame should be returned}
  \item{params}{Three params needed to draw the CI's: fm, distscale and sigline}
  \item{plotPoint}{True or False, indicating whether points without interpolation should be plotted including a CI}
}
\details{
plotPoints defaults to True as there is a certain chance that the area has been trawled / covered / travelled.
}
\references{Hintzen et al. 2010 Improved estimation of trawling tracks using cubic Hermite spline interpolation of position registration data,
EU lot 2 project}
\author{Niels T. Hintzen}
\note{Computation can take a very long time if either grid resolution is high or
if many interpolations are used.
}
\seealso{\code{interpolateVMS(),plotInterpolation()}}
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

  #Plot the interpolated data with a CI
plotCIinterpolation(interpolation,VMS,grid,plot=T,params=list(fm=0.5,distscale=20,sigline=0.2),plotPoint=T)
}
