\name{eflalo2relational}
\alias{eflalo2relational}
\title{
Turn eflalo dataset into relational format
}
\description{
Turn the column setup of eflalo into a row setup where each species catch has its own row
}
\usage{
eflalo2relational(x)
}
\arguments{
  \item{x}{
  Dataframe with eflalo data and eflalo format
}
}
\details{
May take a long time for long eflalo datasets
}

\references{
EU Lot 2 project
}
\author{
Niels T. Hintzen
}

\seealso{
 \code{\link{formatEflalo}},\code{\link{readEflalo}}
}
\examples{
data(eflalo)
eflalo    <- eflalo[1:20,]
eflaloRel <- eflalo2relational(eflalo)
}

