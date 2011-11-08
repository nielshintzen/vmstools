\name{pointInHarbour}
\alias{pointInHarbour}
\title{
Find points in harbour within specified range
}
\description{
Method to find the gps positions given with tacsat data that are situated within a range of a port.
}
\usage{pointInHarbour(lon,lat,harbour)}
\arguments{
  \item{lon}{Longitudinal positions of the TACSAT formatted data}
  \item{lat}{Latitudinal positions of teh TACSAT formatted data}
  \item{harbour}{Latitudinal and Longitudinal position of the harbour and outer range from midpoint of harbour}
}
\details{
The method returns the index of points that are within a harbour area, given the midpoints of the harbours and a 
range (in km) from these midpoints.
}
\references{EU lot 2 project }
\author{Niels T. Hintzen}
\seealso{\code{findEndVMS(),distance(),lonLatRatio(),plotInterpolation(),sortVMS(),
               filterVMS(),mergeEflalo2Tacsat()}}
\examples{
library(vmstools)
data(eflalo)
data(tacsat)
data(euharbours)

#-Remove duplicated records from tacsat
myf       <- paste(tacsat$VE_REF,tacsat$SI_LATI,tacsat$SI_LONG,
                   tacsat$SI_DATE,tacsat$SI_TIME);
tacsat    <- tacsat[!duplicated(myf),];


#-Find all the gps locations that are located within the port area
idx <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,
                      harbours=euharbours))

#-Plot these port locations on a map
library(maps); library(mapdata); library(mapproj);
#map the world, but plot only the northsea by lon and lat limits,
# in high resolution
xrange <- range(tacsat$SI_LONG[idx])
yrange <- range(tacsat$SI_LATI[idx])

map('worldHires',xlim=xrange,ylim=yrange,col="darkgreen",fill=TRUE,
    resolution=1, bg="white", border=0)
map.axes();
mtext("Longitude",1,line=-2,outer=T,cex=1.2,font=2)
mtext("Latitude",2,line=-2,outer=T,cex=1.2,font=2)

points(tacsat$SI_LONG[idx],tacsat$SI_LATI[idx],cex=0.1,pch=16,col="red")
}