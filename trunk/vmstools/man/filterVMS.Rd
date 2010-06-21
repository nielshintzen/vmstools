\name{filterVMS}
\alias{filterVMS}
\title{Filter out duplicated or out-of-range VMS datapoints}
\description{
Filter out duplicated VMS datapoints and possible remove VMS datapoints that are outside
a certain speed range or heading range. A new VMS set, without the out-filtered records is returned.
}
\usage{
filterVMS(VMS, st = c(2, 6), hd = NULL, remDup = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{VMS}{Original VMS data.frame with columns speed, heading, ship and date}
  \item{st}{Speed thresholds in between records are kept. Begin and end threshold or NULL}
  \item{hd}{Heading thresholds in between records are kept. Begin and end threshold or NULL}
  \item{remDup}{Remove duplicates: True or False}
}
\details{
hd is computed as the change in heading from the previous VMS point (if ship remains the same).
Hence, hd does not equal the heading as given in the VMS dataset
}
\references{EU lot 2 project, see also: Mills et al. 2006}
\author{Niels T. Hintzen}
\seealso{\code{sortVMS()}}
\examples{
data(VMS)

  #Sort the VMS data
VMS     <- sortVMS(VMS)
VMS     <- VMS[1:1000,]

  #Filter the VMS data
VMS     <- filterVMS(VMS,c(2,6),hd=NULL,remDup=T) #A warning is given that duplicated VMS points have been removed
}
