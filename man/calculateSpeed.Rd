\name{calculateSpeed}
\alias{calculateSpeed}
\title{
Calculate the speed of the vessel at each ping
}
\description{
Calculate the speed of a vessel based on the time and distance travelled between
pings.
}
\usage{
calculateSpeed(tacsat, level, weight, fill.na)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tacsat}{tacsat dataset}
  \item{level}{level to get calculated speed at: trip or vessel}
  \item{weight}{weight to apply to calculation of mean speed towards and away from ping}
  \item{fill.na}{If speed cannot be calculated based on default or provided
                 weight, take closest alternative to provide a speed}
}
\value{
The original tacsat file is returned including a column: SI_SPCA which holds the calculated speed in knots
}
\details{
Note that the DEFAULT speed given is the difference between ping x and ping x-1. Hence,
the first ping of a vessel or trip does NOT have a calculated speed and will display
NA.

With weight you can specify if the speed is calculated between ping x and ping x-1 (weight = c(1,0)), if
the speed is calculated between ping x and ping x+1 (weight = c(0,1)) or an intermediate weight (weight = c(0.4,0.6) /  equal weight = c(0.5,0.5)).
}

\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{interpolateTacsat}},\code{\link{intervalTacsat}}}
\examples{
data(tacsat)
result <- calculateSpeed(tacsat[1:100,],level="vessel")
result <- calculateSpeed(tacsat[1:100,],level="vessel",weight=c(2,1),fill.na=TRUE)

data(eflalo)
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
result  <- calculateSpeed(tacsatp[1:100,],level="trip",weight=c(1,1),fill.na=FALSE)


}