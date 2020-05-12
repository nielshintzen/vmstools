\name{calculateCI}
\alias{calculateCI}
\title{Calculate the CI between two succeeding VMS datapoints}
\description{
The interpolated tracks can be surrounded by a sort of confidence interval representing
the outer region a vessel could have travelled between two succeeding datapoints. Within
this function the CI's are computed.

}
\usage{    
calculateCI(int,tacint,params,grid,plot=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{int}{interpolation, as data.frame from 'interpolateTacsat'}
  \item{tacint}{tacsat records (two rows) corresponding with the interpolation}
  \item{params}{list of parameters used to perform interpolation}
  \item{grid}{object of class 'GridTopology' specifying the grid dimensions}
  \item{plot}{Logical. Whether the result of the interpolation must be plotted}
}
\value{
Returns the Confidence Interval on a grid of class 'SpatialGridDataFrame' with the
CI values in the data slot.
}
\references{Hintzen et al. 2010 Improved estimation of trawling tracks using cubic Hermite spline interpolation of position registration data,
EU lot 2 project}
\author{Niels T. Hintzen}
\note{Computation can take a very long time if either grid resolution is high or
if many interpolations are used.
}
\seealso{\code{\link{interpolateTacsat}},\code{\link{maxRangeCI}}}
\examples{
data(tacsat)

  #Sort the Tacsat data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Filter the Tacsat data
tacsat          <- filterTacsat(tacsat,c(2,6),hd=NULL,remDup=TRUE)

  #Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
                    res=100,method="cHs",params=list(fm=0.5,distscale=20,
                    sigline=0.2,st=c(2,6)),headingAdjustment=0)

  #Create the final grid where all interpolations should fit on
xrange        <- c(2,3); yrange <- c(51,52)
grid          <- createGrid(xrange,yrange,resx=0.01,resy=0.005)

res           <- calculateCI(interpolation[[4]],
                             tacsat[interpolation[[4]][1,],],
                             params=list(fm=0.25,distscale=3.1,sigline=0.4,st=c(2,6)),
                             grid=grid,
                             plot=TRUE)
}
