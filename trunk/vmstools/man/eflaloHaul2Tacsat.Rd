\name{eflaloHaul2Tacsat}
\alias{eflaloHaul2Tacsat}
\title{
Convert eflalo (haul by haul / log-event) to tacsat format
}
\description{
Convert haul by haul / log event organised eflalo data to tacsat relational style
}
\usage{
eflaloHaul2Tacsat(x, npoints = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
eflalo formatted dataset with haul-by-haul data (e.g. LE_STIME,LE_SLON)
}
  \item{npoints}{
Number of intermediate points to include (defaults to 10)
}
}
\details{
}
\value{
Tacsat formatted dataset is returned
}
\author{
Niels T. Hintzen}
\seealso{
\code{\link{eflalo2Relational}}
}
