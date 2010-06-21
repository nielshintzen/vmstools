\name{sortVMS}
\alias{sortVMS}
\title{Sorting VMS dataframe}
\description{
Sort the VMS data first by vessel, then by date, speed and heading. Needs to be in this order
to be effectively used in other EU lot 2 project generic functions.
}
\usage{
sortVMS(VMS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{VMS}{VMS dataframe}
}

\references{EU lot 2 project}
\author{Niels T. Hintzen}
\note{Uses library(doBy)
}
\seealso{Library: doBY, \code{filterVMS(),data(VMS)}}
\examples{
data(VMS)

  #Sort the VMS data
VMS     <- sortVMS(VMS)
}