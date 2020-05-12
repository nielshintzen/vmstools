\name{pointInHarbour}
\alias{pointInHarbour}
\title{
Find points in harbour within specified range
}
\description{
Method to find the gps positions given with tacsat data that are situated within a range of a port.
}
\usage{pointInHarbour(lon,lat,harbour,rowSize = 30, returnNames = FALSE,
                 saveHarbourList = TRUE)}
\arguments{
  \item{lon}{Longitudinal positions of the TACSAT formatted data}
  \item{lat}{Latitudinal positions of teh TACSAT formatted data}
  \item{harbour}{Latitudinal and Longitudinal position of the harbour and outer range from midpoint of harbour}
  \item{returnNames}{Logical: return the name of the harbour instead of 1 / 0 indicating if it is inside or
  outside the harbour. Default to FALSE}
  \item{saveHarbourList}{Logical: writing harbour list used in function to file. Default to TRUE}
}
\details{
The method returns the index of points that are within a harbour area, given the midpoints of the harbours and a 
range (in km) from these midpoints.
}
\references{EU lot 2 project }
\author{Niels T. Hintzen}
\seealso{\code{\link{distance}}, \code{\link{lonLatRatio}}, \code{\link{sortTacsat}},
         \code{\link{filterTacsat}}, \code{\link{mergeEflalo2Tacsat}}}
\examples{
data(eflalo)
data(tacsat)
data(euharbours); euharbours <- harbours

#-Remove duplicated records from tacsat
myf       <- paste(tacsat$VE_REF,tacsat$SI_LATI,tacsat$SI_LONG,
                   tacsat$SI_DATE,tacsat$SI_TIME);
tacsat    <- tacsat[!duplicated(myf),];


#-Find all the gps locations that are located within the port area
idx <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,
                      harbours=harbours,returnNames=TRUE)
print(head(idx))
getwd() #in this directory, the harbour list will be written to disk
idx <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,
                      harbours=harbours,saveHarbourList=TRUE)

idx <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,
                      harbours=harbours)
idx <- which(idx==1)

#-Plot these port locations on a map
library(maps); library(mapdata)
#map the world, but plot only the northsea by lon and lat limits,
# in high resolution
xrange <- range(tacsat$SI_LONG[idx])
yrange <- range(tacsat$SI_LATI[idx])

map('worldHires',xlim=xrange,ylim=yrange,col="darkgreen",fill=TRUE,
    resolution=1, bg="white", border=0)
map.axes();
mtext("Longitude",1,line=-2,outer=TRUE,cex=1.2,font=2)
mtext("Latitude",2,line=-2,outer=TRUE,cex=1.2,font=2)

points(tacsat$SI_LONG[idx],tacsat$SI_LATI[idx],cex=0.1,pch=16,col="red")
}