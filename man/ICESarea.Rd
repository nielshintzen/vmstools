\name{ICESarea}
\alias{ICESarea}
\title{
Get ICES area from coordinates
}
\description{
Get the ICES area from any lon,lat position, given that this position is
within the ICES region.
}
\usage{
ICESarea(tacsat,areas,proj4string=NULL,fast=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tacsat}{
dataframe given that they have 'SI_LONG' and 'SI_LATI' columns (either tacsat
format or other dataset with SI_LONG and SI_LATI columns)
}
  \item{areas}{
ICES areas as SpatialPolygons}
\item{proj4string}{Projection string, default to NULL.}
\item{fast}{If memory allocation is not a problem, a faster version can be
switched on}
}

\value{
Returns the areas as a vector}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{
\code{\link{ICESrectangle}}, \code{\link{ICESrectangle2LonLat}}
}
\examples{
data(ICESareas)
res   <- data.frame(SI_LONG = c(1,2,2,4,2),
                    SI_LATI = c(53,53.2,54,56.7,55.2))
areas <- ICESarea(res,ICESareas)
}