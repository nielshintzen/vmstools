\name{equalDistance}
\alias{equalDistance}
\title{
Interpolated points at equal distance
}
\description{
Returns the interpolated dataset with only those points remaining that are at equal eucledian
distance from each other, with the number of points to retreive remaining.
}
\usage{equalDistance(interpolation,res=10)}
\arguments{
  \item{interpolation}{interpolated dataset obtained from the interpolation() function}
  \item{res}{number of points to retreive from function}
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{filterTacsat}}, \code{\link{tacsat}}, \code{\link{interpolateTacsat}}}
\examples{
data(tacsat)

#Sort the VMS data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

#Filter the Tacsat data
tacsat     <- filterTacsat(tacsat,st=c(2,6),hd=NULL)

#Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
                    res=100,method="cHs",params=list(fm=0.5,distscale=20,
                    sigline=0.2,st=c(2,6)),headingAdjustment=0)

#Get a set back with only 10 points per interpolation at equal distance
ed_interpolation <- equalDistance(interpolation,10)
}
