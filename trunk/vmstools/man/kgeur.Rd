\name{kgeur}
\alias{kgeur}
\title{
Return kg and euro column index of eflalo dataset
}
\description{
Returns the index of the columns with kg and euro information from a given
eflalo dataset
}
\usage{
kgeur(colnames(eflalo))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Colnames of eflalo dataset
}
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\examples{
data(eflalo)
kgeur(colnames(eflalo))
}

