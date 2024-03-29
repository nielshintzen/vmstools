% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICESrectangle2CSquare.r
\name{ICESrectangle2CSquare}
\alias{ICESrectangle2CSquare}
\title{Convert ICES rectangle to CSquares}
\usage{
ICESrectangle2CSquare(rectangles, degrees, onLand = T)
}
\arguments{
\item{rectangles}{vector with ICES rectangle names}

\item{degrees}{Resolution of CSquare notation: 10, 5, 1, 0.5, 0.1, 0.05,
0.01}

\item{onLand}{logical. Should CSquares with their midpoints on land be returned}
}
\value{
Returns the rectangles as a vector
}
\description{
Give the CSquares belonging to ICES rectangles.
}
\examples{

degrees       <- 0.05
rectangles    <- c("33F3","33F4","32F3","34F3")
squares       <- ICESrectangle2CSquare(rectangles,degrees)
squares       <- do.call(c,squares)
par(xaxs="i",yaxs="i",las=1,oma=c(2,2,1,1))
plot(st_geometry(ICESareas),xlim=c(2,5),ylim=c(51.5,53),col="lightblue",fill=T);
axis(1); axis(2)
abline(v=seq(2,5,1),lty=3)
abline(h=seq(51.5,53,0.5),lty=3)
lonLatSquares <- CSquare2LonLat(squares,degrees=degrees)
points(lonLatSquares$SI_LONG,lonLatSquares$SI_LATI,pch=19,cex=0.2,col=2)
squares       <- ICESrectangle2CSquare(rectangles,degrees,onLand=F)
squares       <- do.call(c,squares)
lonLatSquares <- CSquare2LonLat(squares,degrees=degrees)
points(lonLatSquares$SI_LONG,lonLatSquares$SI_LATI,pch=19,cex=0.2,col=4)


}
\references{
EU Lot 2 project
}
\seealso{
\code{\link{CSquare}}, \code{\link{ICESrectangle2LonLat}}
}
\author{
Niels T. Hintzen
}
