\name{interpolateVMS}
\alias{interpolateVMS}
\title{Find interpolated points between two succeeding GPS datapoints}
\description{
Method to find interpolated points between two succeeding VMS / GPS datapoints. Two different
methods can be applied: the cubic Hermite spline method as developed in Hintzen et al. 2010 or
an interpolation of a straight line (on a x-y plane, not a sphere)
}
\usage{
interpolateVMS(VMS, interval = 120, margin = 12, res = 100, method = "cHs", params = list(fm=0.2,distscale=20,sigline=0.2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{VMS}{VMS dataset with GPS points that need to be interpolated}
  \item{interval}{Time in minutes between the succeeding datapoints}
  \item{Deviation from specified interval to find succeeding datapoints}{}
  \item{res}{Number of points to use to create interpolation (including start and end point), default = 100}
  \item{method}{"cHs" for cubic Hermite spline or "SL" for Straight Line interpolation, default = 0.2,20,0.2}
  \item{params}{Params needed for cHs method}
}
\details{
Interval: In most VMS datasets the succeeding datapoint can be found 1 or 2 hours appart. This interval time
should be specified here. Interval can also be specified as e.g. 15 minutes if the VMS / GPS dataset allows this.
Margin: Hardly ever, the interval time is precise. To allow some deviation from the 1 or 2 hour interval
the margin can be adjusted.

The 'res' value defines the number of points to use to construct the interpolation. Default = 100 which results in a
high resolution interpolation.

Params can be obtained from analyses (see Hintzen et al. 2010) or can be varied to explore the differences. To interpolate
only first value needs to be specified.

The object returned is a list. Each list item contains in the first row the index of the succeeding VMS pionts interpolated.
The second till last row contain the x (first column) and y (second column) interpolated values. The length of the list
indicates the number of interpolations created.
}
\references{Hintzen et al. 2010 Improved estimation of trawling tracks using cubic Hermite spline interpolation
of position registration data, EU lot 2 project }
\author{Niels T. Hintzen}
\seealso{\code{findEndVMS(),distance(),lonLatRatio(),plotInterpolation(),sortVMS(),filterVMS()}}
\examples{
  #Load the VMS data
data(VMS)
  #Filter the VMS data
VMS     <- filterVMS(VMS,st,hd,remDup)

  #Sort the VMS data
VMS     <- sortVMS(VMS)
VMS     <- VMS[1:1000,]
  #Interpolate the VMS data
interpolation <- interpolateVMS(VMS,interval=120,margin=10,res=100,method="cHs",params=list(fm=0.5,distscale=20,sigline=0.2))
}
