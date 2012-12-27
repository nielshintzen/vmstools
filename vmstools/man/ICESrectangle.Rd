\name{ICESrectangle}
\alias{ICESrectangle}
\title{
Get ICES rectangle from coordinates
}
\description{
Get the ICES rectangle from any lon,lat position, given that this position is
within the ICES region.
}
\usage{
ICESrectangle(dF)
}
\arguments{
  \item{dF}{
dataframe given that they have 'SI_LONG' and 'SI_LATI' columns
}
}
\value{
Returns the rectangles as a vector}
\references{EU Lot 2 project}
\author{Neil Campbell}
\seealso{
\code{\link{ICESarea}}, \code{\link{ICESrectangle2LonLat}}}
\examples{
res <- data.frame(SI_LONG = c(1,2,2,4,2),
                  SI_LATI = c(53,53.2,54,56.7,55.2))
ICESrectangle(res)

}