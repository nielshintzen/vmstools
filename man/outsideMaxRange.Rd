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
outsideMaxRange(int,tacint,params,grid)
}
\arguments{
  \item{int}{interpolation, as data.frame from 'interpolateTacsat'}
  \item{tacint}{tacsat records (two rows) corresponding with the interpolation}
  \item{params}{list of parameters used to perform interpolation}
  \item{grid}{object of class 'GridTopology' specifying the grid dimensions}
}
\value{
Returnes a list with three objects:
1) Fraction of CI located outside maximum range
2) Franction of CI located inside maximum range
3) Maximum value of CI (top, at VMS points)

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
tacsat          <- filterTacsat(tacsat,c(2,6),hd=NULL,remDup=TRUE)

  #Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
                    res=100,method="cHs",params=list(fm=0.5,distscale=20,
                    sigline=0.2,st=c(2,6)),headingAdjustment=0)

  #Create the final grid where all interpolations should fit on
xrange        <- c(2,3); yrange <- c(51,52)
grid          <- createGrid(xrange,yrange,resx=0.01,resy=0.005)

res           <- outsideMaxRange(interpolation[[4]],
                                 tacsat[interpolation[[4]][1,],],
                                 params=list(fm=0.25,distscale=3.1,
                                             sigline=0.4,st=c(2,6)),
                                 grid=grid)

}
