\name{pointOnLand}
\alias{pointOnLand}
\title{
Find points on land given a set of coordinates}
\description{
Find the points that are on land given a set of coordinates and polygons that
determine the part that is land
}
\usage{
pointOnLand(tacsat, lands)
}
\arguments{
  \item{tacsat}{
Tacsat file
}
  \item{lands}{
Polygon of area that is considered to be land
}
}
\details{
Function uses the inPoly function. However, with many coordinates, the checking
might take longer.
}
\value{
Returns a vector with values 0 and 1. 1 indicating points on land, 0 indicating points not on land.
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{
\code{pointInHarbour()}
}
\examples{

data(tacsat)
tacsat  <- tacsat[1:1000,]
tacsat  <- sortTacsat(tacsat)

europe          <- list()
for(iLand in c("Netherlands","Germany","Belgium","France","UK",
               "Denmark","Ireland","Norway","Sweden","Iceland")){
  europe[[iLand]] <- map("worldHires",iLand,main=iLand,plot=F,res=0);
}
europe[["Norway"]]$y[which(europe[["Norway"]]$y < 50)]  <- NA
europe[["UK"]]$y    [which(europe[["UK"]]$y < 45)]      <- NA;
  europe[["UK"]]$x    [which(europe[["UK"]]$x < -10)]   <- NA
europe[["France"]]$y[which(europe[["France"]]$y < 40)]  <- NA;
  europe[["France"]]$x[which(abs(europe[["France"]]$x) > 10)]  <- NA

idx             <- pointOnLand(tacsat,europe);
}