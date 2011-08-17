\name{euharbours}
\alias{euharbours}
\docType{data}
\title{
Dataframe containing > 400 harbours primarily Europe, subset of euharbours
}
\description{
A list of harbour positions primarily located in Europe, including name, latitude, longitude and a certain range
the harbour spans
}
\usage{data(harbours)}
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
data(harbours)
pointInHarbour(tacsat$SI_LONG,tacsat$SI_LATI,harbours)
}
\keyword{datasets}