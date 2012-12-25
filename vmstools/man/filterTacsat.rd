\name{filterTacsat}
\alias{filterTacsat}
\title{Filter out duplicated or out-of-range Tacsat datapoints}
\description{
Filter out duplicated Tacsat datapoints and possible remove Tacsat datapoints that are outside
a certain speed range or heading range. A new Tacsat set, without the out-filtered records is returned.
}
\usage{
filterTacsat(tacsat, st = c(2, 6), hd = NULL, remDup = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tacsat}{Original VMS data.frame with columns speed, heading, ship and date}
  \item{st}{Speed thresholds in between records are kept. Begin and end threshold or NULL}
  \item{hd}{Heading thresholds in between records are kept. Begin and end threshold or NULL}
  \item{remDup}{Remove duplicates: True or False}
}
\details{
hd is computed as the change in heading from the previous Tacsat point (if ship remains the same).
Hence, hd does not equal the heading as given in the Tacsat dataset
}
\references{EU lot 2 project, see also: Mills et al. 2006}
\author{Niels T. Hintzen}
\seealso{\code{\link{sortTacsat}}}
\examples{
data(tacsat)

  #Sort the VMS data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Filter the VMS data
  #A warning is given that duplicated VMS points have been removed
tacsat     <- filterTacsat(tacsat,c(2,6),hd=NULL,remDup=TRUE)
}
