\name{effort}
\alias{effort}
\title{
Calculate effort from tacsat or eflalo file
}
\description{
Calculate the effort in the tacsat or eflalo file by trip
}
\usage{
effort(x, level = "trip", unit = "hours",weight=c(0.5,0.5),fill.na=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
tacsat or eflalo object
}
  \item{level}{
"trip"
}
  \item{unit}{
units in c("auto", "secs", "mins", "hours",
                   "days", "weeks")
}
  \item{weight}{weight to apply to calculation of mean speed towards and away from ping}
  \item{fill.na}{If speed cannot be calculated based on default or provided
                 weight, take closest alternative to provide a speed}

}
\details{
Note that the DEFAULT effort given is based on the difference between ping x and
ping x-1. Hence, the first ping of a vessel or trip does NOT have a calculated
speed and will display NA.

With weight you can specify if the speed is calculated between ping x and ping x-1
(weight = c(1,0)), if the speed is calculated between ping x and ping x+1 (weight = c(0,1))
or an intermediate weight (weight = c(0.4,0.6) /  equal weight = c(0.5,0.5)).
}
\value{
Returnes the original tacsat or eflalo file but including an effort column
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\examples{
data(tacsat)
data(eflalo)

tacsat            <- sortTacsat(tacsat)
tacsatp           <- mergeEflalo2Tacsat(eflalo,tacsat)
teffort           <- effort(tacsatp,level="trip",unit="hours")
eeffort           <- effort(eflalo,level="trip",unit="hours")
}