\name{euharbours}
\alias{euharbours}
\alias{harbours}
\docType{data}
\title{
Dataframe containing > 3500 harbours primarily Europe
}
\description{
A list of harbour positions primarily located in Europe, including
-name,
-latitude,
-longitude and a certain
-range (in km) the harbour spans (radius).
}
\usage{data(euharbours)}
\references{
EU Lot 2 project, Niels T. Hintzen (niels.hintzen@wur.nl)
}
\examples{
data(tacsat)
data(euharbours); euharbours <- harbours

#modify the range of the harbours to e.g. 2km
euharbours$range <- 2

pointInHarbour(tacsat$SI_LONG,tacsat$SI_LATI,euharbours)
}
\keyword{datasets}
