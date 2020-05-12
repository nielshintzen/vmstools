\name{an}
\alias{an}
\title{shortcut for as.numeric}
\description{
Change the class of an object to numeric
}
\usage{
an(x)
}
\arguments{
  \item{x}{object to turn into numeric
}}
\value{
as.numeric attempts to coerce its argument to numeric type
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{ \code{\link{as.numeric}}}
\examples{
as.numeric("5")   #returns the character 5 as class 'numeric'
an("5")           #returns the character 5 also as class 'numeric'
}
