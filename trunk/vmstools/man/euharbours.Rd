\name{euharbours}
\alias{euharbours}
\docType{data}
\title{
Dataframe containing > 3500 harbours primarily Europe
}
\description{
A list of harbour positions primarily located in Europe, including name, latitude, longitude and a certain range
the harbour spans
}
\usage{data(euharbours)}
\format{
}
\details{
}
\source{
}
\references{
EU Lot 2 project, Niels T. Hintzen (niels.hintzen@wur.nl)
}
\examples{
data(euharbours)
pointInHarbour(tacsat$SI_LONG,tacsat$SI_LATI,euharbours)
}
\keyword{datasets}
