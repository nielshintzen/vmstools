\name{sortTacsat}
\alias{sortTacsat}
\title{Sorting Tacsat dataframe}
\description{
Sort the Tacsat data first by vessel, then by date, speed and heading. Needs to be in this order
to be effectively used in other EU lot 2 project generic functions.
}
\usage{
sortTacsat(tacsat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tacsat}{tacsat dataframe}
}

\references{EU lot 2 project}
\author{Niels T. Hintzen}
\note{Uses library(doBy)
}
\seealso{Library: doBY, \code{filterTacsat(),data(tacsat)}}
\examples{
data(tacsat)

  #Sort the VMS data
tacsat     <- sortTacsat(tacsat)
}