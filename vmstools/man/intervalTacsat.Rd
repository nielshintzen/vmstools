\name{intervalTacsat}
\alias{intervalTacsat}
\title{
Return the interval time between pings
}
\description{
Return the interval time between pings of one vessel or at the trip level. Interval
is calculated based on ping x and ping x-1 or ping x and ping x+1.
}
\usage{
intervalTacsat(tacsat, level,weight,fill.na)
}
\arguments{
  \item{tacsat}{tacsat dataset}
  \item{level}{level to get interval rate at: trip or vessel}
  \item{weight}{weight to apply to calculation of mean interval rate towards and away from ping}
  \item{fill.na}{If interval rate cannot be calculated based on default or provided
                 weight, take closest alternative to provide an interval rate}
}
\value{
The original tacsat file is returned including a column: INTV which holds the interval rate in minutes
}
\details{
Note that the DEFAULT interval given is the difference between ping x and ping x-1. Hence,
the first ping of a vessel or trip does NOT have an interval rate and will display
NA.

With weight you can specify if the interval rate is used between ping x and ping x-1 (weight = c(1,0)), if
the interval rate is used between ping x and ping x+1 (weight = c(0,1)) or an intermediate weight (weight = c(0.4,0.6) /  equal weight = c(0.5,0.5)).
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{interpolateTacsat}},\code{\link{calculateSpeed}}}
\examples{
data(tacsat)
result <- intervalTacsat(tacsat[1:100,],level="vessel")
result <- intervalTacsat(tacsat[1:100,],level="vessel",weight=c(2,1),fill.na=TRUE)

data(eflalo)
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
result  <- intervalTacsat(tacsatp[1:100,],level="trip",weight=c(1,1),fill.na=FALSE)

}

