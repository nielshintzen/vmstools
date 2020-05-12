\name{polyDef}
\alias{polyDef}
\title{
Define polygon in text string format}
\description{
Create a text string format polygon which can be called by other functions}
\usage{
polyDef(lon, lat, gridx, gridy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lon}{
Vector with longitude (or x) values
}
  \item{lat}{
Vector with latitude (or y) values, should have the same length as \code{lon}
}
  \item{gridx}{
gridcell size in degrees in the longitude / x direction
}

  \item{gridy}{
gridcell size in degrees in the latitude / y direction
}
}
\details{
The function stets up a grid - with the origin at (0,0) - and assigns a grid cell to each of the points given by \code{lon} and \code{lat}.
}
\value{
A Well Known Text string for each value of \code{lon} and \code{lat}.
}
\note{
Any points that lie exactly on the border between two grid cells will be assigned to the grid cell above in the northern hemisphere, below in the southern hemisphere, to the right in the eastern hemisphere and to the left in the western hemisphere.
}
\author{
Hans D Gerritsen
}
\seealso{
\code{\link{nestedGrid}}}
\examples{
lon <- rnorm(25,3)
lat <- rnorm(25,53)
pols <- polyDef(lon, lat, gridx = 0.5, gridy= 0.5)
plot(lon,lat)
tempfun <- function(i){polygon(eval(parse(text=pols[i]))@coords)}
lapply(as.list(1:length(pols)),tempfun)
}
