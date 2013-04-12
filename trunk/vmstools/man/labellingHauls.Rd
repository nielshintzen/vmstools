\name{labellingHauls}
\alias{labellingHauls}
\title{Labelling fishing operations}
\description{
Labelling fishing sequences from the SI_STATE coding with an unique identifier
}
\usage{
labellingHauls(tacsat)
}
\arguments{
  \item{tacsat}{ a tacsat format with the SI_STATE coding informed for "h", "f" and "s"
}
}
\value{
add a new column to the tacsat data.frame named HL_ID
}
\references{}
\author{F. Bastardie}
\seealso{\code{\link{activityTacsat}}}
\examples{
data(tacsat)
tacsat$SI_STATE <- 0
tacsat$SI_SP    <- replace(tacsat$SI_SP, is.na(tacsat$SI_SP), 0)
tacsat[tacsat$SI_SP >= 1.5 & tacsat$SI_SP <= 7.5,'SI_STATE']      <- 'f'
tacsat[tacsat$SI_SP < 1.5,'SI_STATE']                              <- 'h'
tacsat[tacsat$SI_SP > 7.5,'SI_STATE']                              <- 's'
labellingHauls(tacsat)
}
