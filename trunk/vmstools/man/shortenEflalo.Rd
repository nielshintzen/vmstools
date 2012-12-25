\name{shortenEflalo}
\alias{shortenEflalo}
\title{Cuts up eflalo according to any particular combination of species to make it more manageable
}
\description{
Cuts up eflalo according to any particular combination of species to make it more manageable
}
\usage{
shortenEflalo(data = eflalo, which.species = c("SPEC"))
}
\arguments{
  \item{data}{
eflalo formatted data
}
  \item{which.species}{
array of species names in FAO species names
}
}
\references{
EU lot 2 project
}
\author{
Doug Beare
}
\examples{
data(eflalo)
shortenEflalo(data = eflalo, which.species = c("PLE"))

}