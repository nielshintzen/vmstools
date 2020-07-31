\name{formatTacsat}
\alias{formatTacsat}
\title{Format tacsat data to ensure that all columns are in right format
}
\description{
Reformat all the columns of the tacsat data to ensure that all data is in the right format
}
\usage{
formatTacsat(x)
}
\arguments{
  \item{x}{ tacsat file
}
}
\value{
Returns original tacsat file but now with reformatted data}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\seealso{\code{\link{formatEflalo}}
}
\examples{
data(tacsat)
tacsat <- formatTacsat(tacsat)
}
