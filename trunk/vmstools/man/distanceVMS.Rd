\name{distanceVMS}
\alias{distanceVMS}
\title{Calculate the distance (in km) from a VMS dataset}
\description{
Calculates the distance of a VMS dataset with specification of succeeding datapoints. Distance is only
calculated between these specified datapoints}
\usage{
distanceVMS(VMS, index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{VMS}{VMS (normal or high ping rate) dataset}
  \item{index}{Matrix with specification of succeeding datapoints (see details for format)}
}
\details{
index is designed as a matrix where rows represent succeeding datapoints, column 1 represent start points, column 2 represent end points.
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{distance(),distanceInterpolation(),diffInter()}}
\examples{

data(VMS)
st        <- c(2,6)       #Speed threshold points (two values), NULL means use all points
remDup    <- T            #Remove duplicate records in VMS dataset

  #Sort the VMS data
VMS     <- sortVMS(VMS)
VMS     <- VMS[1:1000,]

  #Filter the VMS data
VMS     <- filterVMS(VMS,st,NULL,remDup)

distanceVMS(VMS,matrix(c(2,3,3,4),nrow=2,ncol=2,dimnames=list(1:2,c("startpoint","endpoint")))) #0.1044712 0.2875203
}

