\name{segmentTacsatSpeed}
\alias{segmentTacsatSpeed}
\title{
do the distinction between fishing/non-fishing detecting the fishing speed peak 
}
\description{
do the distinction between fishing/non-fishing based on (semi)automatic detection 
of the fishing speed peak. This is done per vessel and independantly for each gear
of the vessel. A segmented regression is performed on the cumulative speed distribution
to identify the speed breakpoints bounding the fishing peak.
}
\usage{segmentTacsatSpeed(tacsat,vessels,general)}
\arguments{
  \item{tacsat}{tacsat with LE_GEAR informed for each ping}
  \item{vessels}{a vector of vessel identifiers VE_REF found in tacsat}
  \item{general}{general settings e.g. the output path for saving the plot}
}
\details{
The method returns tacsat with the SI_STATE informed.
}
\references{Bastardie et al. 2010. Fisheries Research}
\author{Francois Bastardie}
\seealso{\code{}}
\examples{
data(tacsat)
tacsat <- segmentSpeedTacsat (tacsat, vessels="35",
              general=list(output.path=file.path('C:','output')))


}