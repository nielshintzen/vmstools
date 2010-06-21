\name{distanceInterpolation}
\alias{distanceInterpolation}

\title{Calculate the distance (in km) from an interpolated dataset}
\description{
For each interpolation between two succeeding datapoints, the distance travelled is computed
based on the function 'distance()'. This function is an easy to use wrapper for a whole
interpolated dataset
}
\usage{
distanceInterpolation(interpolation)
}
\arguments{
  \item{interpolation}{Interpolated dataset as output from the function 'interpolateVMS'}
}
\details{On default, each interpolation consists of 100 points. Hence, 99 distances are computed 
for each interpolation. Returned are 5 }
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{distance(), distanceVMS()}}
\examples{
data(VMS)
st        <- c(2,6)       #Speed threshold points (two values), NULL means use all points
remDup    <- T            #Remove duplicate records in VMS dataset

  #Sort the VMS data
VMS     <- sortVMS(VMS)
VMS     <- VMS[1:1000,]

  #Filter the VMS data
VMS     <- filterVMS(VMS,st,NULL,remDup)

interpolation <- interpolateVMS(VMS,interval=120,margin=10,res=100,method="cHs",params=list(fm=0.2,distscale=20,sigline=0.2))
  #Number of interpolations:
length(interpolation)
  #Calculate distance:
distanceInterpolation(interpolation) #Returns 159 values, which represents the distance travelled in each of the 159 interpolations
}
