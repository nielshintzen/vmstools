\name{sortTacsat}
\alias{sortTacsat}
\title{Sorting Tacsat dataframe}
\description{
Sort the Tacsat data first by vessel, then by date, speed and heading. Needs to be in this order
to be effectively used in other EU lot 2 project generic functions.
}
\usage{
sortTacsat(dat)
}
\arguments{
  \item{dat}{tacsat dataframe}
}

\references{EU lot 2 project}
\author{Niels T. Hintzen}
\note{Uses library(doBy)
}
\seealso{\code{\link{filterTacsat}}}
\examples{
data(tacsat)
require(doBy)

  #Sort the Tacsat data
tacsat     <- sortTacsat(tacsat)
}