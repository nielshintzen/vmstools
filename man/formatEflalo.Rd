\name{formatEflalo}
\alias{formatEflalo}
\title{Format Eflalo data to ensure that all columns are in right format
}
\description{
Reformat all the columns of the Eflalo data to ensure that all data is in the right format
}
\usage{
formatEflalo(x)
}
\arguments{
  \item{x}{ eflalo file
}
}
\value{
Returns original Eflalo file but now with reformatted data}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{formatTacsat}}
}
\examples{
data(eflalo)
eflalo <- formatEflalo(eflalo)
}