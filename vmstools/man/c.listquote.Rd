\name{c.listquote}
\alias{c.listquote}
\title{
Produce quoted list for use within data.table
}
\description{
Produce the quoted list format from a character vector which is used as input for the
data.table format
}
\usage{
c.listquote(...)
}
\arguments{
  \item{\dots}{
character vector
}
}
\details{
character vector could contain e.g. colnames(x) or c("COL1","COL2","COL3") where COL1 etc are column names
}
\value{
Returnes quoted list
}
\references{
See \code{\link{data.table}} documentation
}
\author{
Niels T. Hintzen
}
\note{
With great thanks to original author.
}
\seealso{
See \code{\link{data.table}} documentation
}
\examples{
c.listquote(c("COL1","COL2","COL3")) #Returns: list(COL1,COL2,COL3)
}
