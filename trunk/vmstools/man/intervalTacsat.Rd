\name{intervalTacsat}
\alias{intervalTacsat}
\title{
Return the interval time between pings
}
\description{
Return the interval time between pings of one vessel or at the trip level.
}
\usage{
intervalTacsat(tacsat, level = "trip")
}
\arguments{
  \item{tacsat}{tacsat dataset}
  \item{level}{level to get interval rate at: trip or vessel}
}
\value{
The original tacsat file is returned including a column: INTV which holds the interval rate in minutes
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{interpolateTacsat()}}
\examples{
data(tacsat)

result <- intervalTacsat(tacsat[1:100,],level="vessel")
}

