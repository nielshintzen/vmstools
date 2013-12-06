\name{ICESrectangle2LonLat}
\alias{ICESrectangle2LonLat}
\title{
Convert the ICES rectangle to coordinates
}
\description{
Convert the ICES rectangle to the longitude and latitude position of that
specific rectangle, with either midpoint or lower-left corner.
}
\usage{
ICESrectangle2LonLat(statsq, midpoint = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{statsq}{
character vector of ICES rectangles}
  \item{midpoint}{
logical statement if returned values need to be midpoints of rectangle or not (default)
}
}
\value{
Returns dataframe with longitude and latitude columns}
\references{EU Lot 2 project}
\author{Neil Campbell, Mark Payne}
\seealso{
\code{\link{ICESarea}}, \code{\link{ICESrectangle}}}
\examples{
res       <- data.frame(SI_LONG = c(1,2,2,4,2),
                        SI_LATI = c(53,53.2,54,56.7,55.2))
rects     <- ICESrectangle(res)
midpoints <- ICESrectangle2LonLat(rects,midpoint=TRUE)
corners   <- ICESrectangle2LonLat(rects,midpoint=FALSE)

}