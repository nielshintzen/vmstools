\name{eflalo2Pings}
\alias{eflalo2Pings}
\title{
Merge eflalo to tacsat pings
}
\description{
Internal function of splitAmongPings to merge eflalo landings or values to
tacsat pings
}
\usage{
eflalo2Pings(eflalo, tacsat, pings, vars, eflaloCol, remainTacsat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eflalo}{
eflalo dataset
}
  \item{tacsat}{
tacsat dataset
}
  \item{pings}{
number of pings by variable
}
  \item{vars}{
variable to merge eflalo to tacsat
}
  \item{eflaloCol}{
column names of eflalo
}
  \item{remainTacsat}{
number of tacsat pings that have not been merged yet
}
}
\value{
Returns a list of the eflalo dataset, but without the landings and values that have
been merged, returns the merged tacsat dataset and returns the number tacsat pings that
have not been merged yet}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{
\code{\link{splitAmongPings}}, \code{\link{mergeEflalo2Tacsat}}, \code{\link{mergeEflalo2Pings}}
}
