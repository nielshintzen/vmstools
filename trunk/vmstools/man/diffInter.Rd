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
\seealso{\code{\link{distance}}, \code{\link{distanceInterpolation}}, \code{\link{distanceVMS}}}
\examples{
\dontrun{
data(VMShf)

#-Put the data in the right format
VMShf$SI_DATE <- format(as.Date(VMShf$date),"%d/%m/%Y")
VMShf$SI_TIME <- format(VMShf$date, "%H:%M")

colnames(VMShf) <- c("VE_REF","SI_LATI","SI_LONG","SI_SP",
                     "SI_HE","SI_DATIM","SI_DATE","SI_TIME")

#-Sort the data and remove non-fishing pings
VMShf   <- sortTacsat(VMShf)
VMShf   <- filterTacsat(VMShf,c(2,6),NULL,T)

interpolation <- interpolateTacsat(VMShf,interval=120,margin=10,res=100,
                    method="cHs",params=list(fm=0.3,distscale=20,
                    sigline=0.2,st=c(2,6)),headingAdjustment=0)

#Returns a matrix with 21 rows (each row represents 1 interpolation)
# and 5 measures
cHs <- diffInter(interpolation,VMShf)
}
}
