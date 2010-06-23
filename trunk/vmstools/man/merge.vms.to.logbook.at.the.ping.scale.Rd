\name{merge.vms.to.logbook.at.the.ping.scale}
\alias{merge.vms.to.logbook.at.the.ping.scale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coupling of VMS and logbooks data including landings dispatching between fishing pings.
}
\description{
The first step of the merging routine consisted of assigning a common identifier joining each
VMS trip to a logbook trip by searching for the nearest trip among all logbook trips based on the
temporal mid-point of trips. This circumvents problems where the start and end time of VMS
and logbook trips do not match exactly and ensure that each logbook trip will be merged with a
given VMS trip even if the reporting of the start and end date by fishermen is uncertain. We
observed in the data that this reporting date may fluctuate with +/- 1 day around the VMS trip
dates possibly because of error in reporting. Additionally, the possible few remaining logbook
trips were also linked with the nearest VMS trips. This latter event may occur from unknown
failure in the VMS device while fishermen were still fishing. In the processing, links were then
forced to nearest VMS trips to prevent loss of any reported landings. The second step of the
process evaluated the consistency between both types of data sources. Sequentially, a quality
flag was issued for each recorded position depending on the various degrees of matching, from
finer to coarser, if both data have in common: (i) the trip identifier, the same area (ICES
rectangle), and the catch date (ii) only the trip identifier and the area, or (iii) only the trip
identifier. The initial datasets (both VMS and logbooks) were split into sub-blocks according to
this flag. Each pair of sub-blocks was then joined separately and all merged sub-blocks were
bound in one dataset afterwards.
}
\usage{
merge.vms.to.logbook.at.the.ping.scale(logbooks, general = general, ...)
}

\arguments{
  \item{logbooks}{
%%     eflalo format
}
  \item{general}{
}
  \item{\dots}{
}
}
\details{
in tacsat SI_STAT, fishing should be coded 1 and non-fishing should be coded 2
}
\value{
a merged data.frame is returned
}
\references{
Bastardie et al. 2010
}
\author{
Francois Bastardie
}
\note{
}


\seealso{
}

\examples{

  an <<- function(x) as.numeric(as.character(x)) # alias

data(eflalo)

 merge.vms.to.logbook.at.the.ping.scale (logbooks=eflalo, general=general, 
                 a.vesselid=c("vessel1","vessel2"))
  #=> per vessel, merge logbook with vms

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
