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
ICESarea(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
dataframe given that they have 'SI_LONG' and 'SI_LATI' columns
}
}

\value{
Returns the areas as a vector}
\references{EU Lot 2 project}
\author{Neil Campbell}
\seealso{
\code{ICESrectangle(), ICESrectangle2LonLat()
}
\examples{
res   <- data.frame(SI_LONG = c(1,2,2,4,2),SI_LATI = c(53,53.2,54,56.7,55.2))
areas <- ICESarea(res)
}