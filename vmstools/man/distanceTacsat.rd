\name{distanceTacsat}
\alias{distanceTacsat}
\title{Calculate the distance (in km) from a Tacsat dataset}
\description{
Calculates the distance of a Tacsat dataset with specification of succeeding datapoints. Distance is only
calculated between these specified datapoints}
\usage{
distanceTacsat(tacsat, index)
}
\arguments{
  \item{tacsat}{tacsat (normal or high ping rate) dataset}
  \item{index}{Matrix with specification of succeeding datapoints (see details for format)}
}
\details{
index is designed as a matrix where rows represent succeeding datapoints, column 1 represent start points, column 2 represent end points.
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{distance}}, \code{\link{distanceInterpolation}}, \code{\link{diffInter}}}
\examples{

data(tacsat)
#Speed threshold points (two values), NULL means use all points
st        <- c(2,6)
#Remove duplicate records in VMS dataset
remDup    <- TRUE

  #Sort the VMS data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Filter the VMS data
tacsat     <- filterTacsat(tacsat,st,NULL,remDup)

distanceTacsat(tacsat,matrix(c(2,3,3,4),nrow=2,ncol=2,
                dimnames=list(1:2,c("startpoint","endpoint"))))
#6.335944 14.847291
}

