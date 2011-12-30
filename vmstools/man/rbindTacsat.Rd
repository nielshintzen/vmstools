\name{rbindTacsat}
\alias{rbindTacsat}
\title{
Row bind two tacsat sets together
}
\description{
Row binds two tacsat sets together while taking differences in column names into
account
}
\usage{
rbindTacsat(set1, set2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{set1}{
  Tacsat set 1
}
  \item{set2}{
  Tacsat set 2
}
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{rbindEflalo}}, \code{\link{do.call}}}

\examples{
data(tacsat)
set1 <- tacsat
set2 <- tacsat[seq(1,100,5),]

combined <- rbindTacsat(set1,set2)

}
