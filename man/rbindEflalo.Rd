\name{rbindEflalo}
\alias{rbindEflalo}
\title{
Row bind two eflalo sets together
}
\description{
Row binds two eflalo sets together while taking differences in column names into
account
}
\usage{
rbindEflalo(set1, set2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{set1}{
  Eflalo set 1
}
  \item{set2}{
  Eflalo set 2
}
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{rbindTacsat}}, \code{\link{do.call}}}

\examples{
data(eflalo)
set1 <- eflalo
set2 <- eflalo[seq(1,100,5),]

combined <- rbindEflalo(set1,set2)

}
