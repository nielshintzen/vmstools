\name{ac}
\alias{ac}
\title{shortcut for as.character}
\description{
Change the class of an object to character
}
\usage{
ac(x)
}
\arguments{
  \item{x}{object to turn into character
}
}
\value{
as.character attempts to coerce its argument to character type
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{as.character}}}
\examples{
as.character(5) #returns the number 5 as class 'character'
ac(5)           #returns the number 5 also as class 'character'
}
