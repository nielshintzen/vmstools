\name{af}
\alias{af}
\title{shortcut for as.factor}
\description{
Change the class of an object to factor
}
\usage{
af(x)
}
\arguments{
  \item{x}{object to turn into factor
}}

\value{
as.factor attempts to coerce its argument to factor type
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{ \code{\link{as.factor}}}
\examples{
as.factor(5)    #returns the number 5 as class 'factor'
af(5)           #returns the number 5 also as class 'factor'
}
