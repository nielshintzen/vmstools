\name{interpolateTacsat}
\alias{interpolateTacsat}
\title{Find interpolated points between two succeeding GPS datapoints}
\description{
Method to find interpolated points between two succeeding Tacsat / GPS datapoints. Two different
methods can be applied: the cubic Hermite spline method as developed in Hintzen et al. 2010 or
an interpolation of a straight line (on a x-y plane, not a sphere)
}
\usage{
interpolateTacsat(tacsat, interval = 120, margin = 12, res = 100,
                  method = "cHs", params = list(fm=0.5,distscale=20,
                  sigline=0.2,st=c(2,6)),headingAdjustment=0,fast=FALSE)
}
\arguments{
  \item{tacsat}{Tacsat dataset with GPS points that need to be interpolated}
  \item{interval}{Time in minutes between the succeeding datapoints}
  \item{margin}{Deviation from specified interval to find succeeding datapoints}
  \item{res}{Number of points to use to create interpolation (including start and end point), default = 100}
  \item{method}{"cHs" for cubic Hermite spline or "SL" for Straight Line interpolation, default = 0.2,20,0.2}
  \item{params}{Params needed for cHs method, including the speed thresholds used}
  \item{headingAdjustment}{Parameter to adjust the choice of heading depending on its own or previous point}
  \item{fast}{If your machine holds enough RAM (4gig or more) than the fast version
              runs 7x faster}
}
\details{
Interval: In most Tacsat datasets the succeeding datapoint can be found 1 or 2 hours appart. This interval time
should be specified here. Interval can also be specified as e.g. 15 minutes if the Tacsat / GPS dataset allows this.
Margin: Hardly ever, the interval time is precise. To allow some deviation from the 1 or 2 hour interval
the margin can be adjusted.

The 'res' value defines the number of points to use to construct the interpolation. Default = 100 which results in a
high resolution interpolation.

Params can be obtained from analyses (see Hintzen et al. 2010) or can be varied to explore the differences. To interpolate
only first value needs to be specified and the speedtreshold needs to be given.

headingAdjustment can either be 0 or 1 in case the heading at the endpoint does not represent the heading of the arriving vessel to that
point but the departing vessel. This is only of use when heading needs to be calculated outside the GPS dataset where heading is calculated
at the bearing between two points. For the startpoint, the bearing is correct, for the endpoint however, the heading is incorrect. 

The object returned is a list. Each list item contains in the first row the index of the succeeding Tacsat pionts interpolated.
The second till last row contain the x (first column) and y (second column) interpolated values. The length of the list
indicates the number of interpolations created.
}
\references{Hintzen et al. 2010 Improved estimation of trawling tracks using cubic Hermite spline interpolation
of position registration data, EU lot 2 project }
\author{Niels T. Hintzen}
\seealso{\code{\link{findEndTacsat}},
         \code{\link{distance}},
         \code{\link{lonLatRatio}},
         \code{\link{plotInterpolation}},
         \code{\link{sortTacsat}},
         \code{\link{filterTacsat}}}
\examples{
  #Load the Tacsat data
data(tacsat)

  #Sort the VMS data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Filter the Tacsat data
tacsat     <- filterTacsat(tacsat,st=c(2,6),hd=NULL)

  #Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,res=100,
                                   method="cHs",params=list(fm=0.5,distscale=20,
                                   sigline=0.2,st=c(2,6)),headingAdjustment=0)
}
