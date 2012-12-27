\name{CSquare2LonLat}
\alias{CSquare2LonLat}
\title{Calculate the GPS notation from the CSquare notation

}
\description{
Convert the CSquare notation to GPS location where you can define what resolution
you want the GPS notation to be.
}
\usage{
CSquare2LonLat(csqr, degrees)
}
\arguments{
  \item{csqr}{
CSquare notation like: 1400:134:120:356
}
  \item{degrees}{
Degrees resolution of returned GPS notation : 10, 5, 1, 0.5, 0.1, 0.05, 0.01
}
}
\references{EU Lot 2 project, based on CSquare: http://www.marine.csiro.au/csquares/ }
\author{Niels T. Hintzen}
\seealso{\code{\link{km2Degree}}, \code{\link{degree2Km}}, \code{\link{lonLatRatio}}, \code{\link{distance}}}
\examples{
data(tacsat)
tacsat$CSquare  <- CSquare(tacsat$SI_LONG,tacsat$SI_LATI,0.05)
tacsat$SI_LONG2 <- CSquare2LonLat(tacsat$CSquare,0.5)$SI_LONG
tacsat$SI_LATI2 <- CSquare2LonLat(tacsat$CSquare,0.5)$SI_LATI
tacsat$CSquare2 <- CSquare(tacsat$SI_LONG2,tacsat$SI_LATI2,1)
head(tacsat[which(tacsat$CSquare != tacsat$CSquare2),])
}
