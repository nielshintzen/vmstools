\name{diffInter}
\alias{diffInter}
\title{Calculate distance between interpolation and reference set}
\description{Calculate for each interpolation the distance at fixed points, depending on
the number of points present in the reference set, between the interpolated trajectory
and the reference trajectory. This indicates the deviation of the interpolated set to the
reference set. 
}
\usage{
diffInter(interpolation, reference)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interpolation}{Interpolated dataset as output from the function 'interpolateVMS'}
  \item{reference}{Original high higher resolution VMS dataset}
}
\details{Each interpolation has a start and end point which are similar to the start and end point
of the reference set. In between these two points, the reference set can have more in-between points which are
not used for interpolation, depending on the interval chosen in the 'interpolateVMS' function. These points are matched up with
the points on the interpolated track. For each of these reference in-between points, it is calculated which fraction of the total
distance was travelled. The same fraction is applied to the interpolated track, and the point that matches with this distance is 
linked to the reference set point. The distance function is used to calculate the distance between the reference and interpolated 
dataset. 

The mean, logmean, standard deviation, sdlog and total sum of these distances is returned as a matrix.
}
\references{Hintzen et al. 2010 Improved estimation of trawling tracks using cubic Hermite spline interpolation of position registration data,
EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{distance(),distanceInterpolation(),distanceVMS()}}
\examples{
data(VMShf)

VMShf   <- sortVMS(VMShf)
VMShf   <- filterVMS(VMShf,c(2,6),NULL,T)

interpolation <- interpolateVMS(VMShf,interval=120,margin=5,res=100,method="cHs",params=list(fm=0.2,distscale=20,sigline=0.2))

cHs <- diffInter(interpolation,VMShf) #Returns a matrix with 38 rows (each row represents 1 interpolation) and 5 measures
}
