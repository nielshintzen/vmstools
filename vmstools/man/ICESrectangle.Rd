% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICESrectangle.r
\name{ICESrectangle}
\alias{ICESrectangle}
\title{Get ICES rectangle from coordinates}
\usage{
ICESrectangle(dF)
}
\arguments{
\item{dF}{dataframe given that they have 'SI_LONG' and 'SI_LATI' columns}
}
\value{
Returns the rectangles as a vector
}
\description{
Get the ICES rectangle from any lon,lat position, given that this position
is within the ICES region.
}
\examples{

res <- data.frame(SI_LONG = c(1,2,2,4,2),
                  SI_LATI = c(53,53.2,54,56.7,55.2))
ICESrectangle(res)


}
\references{
EU Lot 2 project
}
\seealso{
\code{\link{ICESarea}}, \code{\link{ICESrectangle2LonLat}}
}
\author{
Neil Campbell
}
