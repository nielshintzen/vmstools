\name{effort}
\alias{effort}
\title{
Calculate effort from tacsat or eflalo file
}
\description{
Calculate the effort in the tacsat or eflalo file by trip
}
\usage{
effort(x, level = "trip", unit = "hours")
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
}
\details{
Currently only level="trip" is implemented
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