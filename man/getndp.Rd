\name{getndp}
\alias{getndp}
\title{Get Number of Decimal Places}
\description{
Return the number of decimal places of a 'numeric'.
}
\usage{
getndp(x, tol = 2 * .Machine$double.eps)
}
\arguments{
  \item{x}{Number to find decimal places from}
  \item{tol}{Tolerance to use}
}
\references{EU lot 2 project}
\author{See R-help pages}
\note{Function not created under EU lot 2 project but found under R-help. Please look there
for credits.
}
\examples{
getndp(5.677) #result: 3
}
