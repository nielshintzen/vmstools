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
\seealso{\code{\link{distance}}, \code{\link{distanceTacsat}}}
\examples{
data(tacsat)
#Speed threshold points (two values), NULL means use all points
st        <- c(2,6)
#Remove duplicate records in VMS dataset
remDup    <- TRUE

  #Sort the tacsat data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Filter the tacsat data
tacsat     <- filterTacsat(tacsat,st,NULL,remDup)

interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
                  res=100,method="cHs",params=list(fm=0.3,distscale=20,
                  sigline=0.2,st=c(2,6)),headingAdjustment=0)
  #Number of interpolations:
length(interpolation)
  #Calculate distance:
  #Returns 159 values, which represents the distance travelled in
  #each of the 159 interpolations
distanceInterpolation(interpolation)
}
