\name{tacsatMinimumConvexPolygon}
\alias{tacsatMinimumConvexPolygon}
\title{
Flag tacsat records that are within the convex polygon
}
\description{
Flag tacsat records that are within the convex polygon with a predefined threshold
}
\usage{
tacsatMinimumConvexPolygon(tacsat, pctThreshold = 90)
}
\arguments{
  \item{tacsat}{
  Tacsat dataframe
}
  \item{pctThreshold}{
  Threshold of points to consider. Between 0 and 100.
}
}
\details{
See point.in.polygon function for returned value details}
\references{EU Lot 2 project}
\author{
Fabrizio Manco
}

\examples{
require(adehabitat)
require(ade4)
data(tacsat)
tacsat <- tacsat[1:100,]
tacsatMinimumConvexPolygon(tacsat,95)
}
