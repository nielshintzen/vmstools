\name{segmentTacsatSpeed}
\alias{segmentTacsatSpeed}
\title{         
do the distinction between fishing/non-fishing detecting the fishing speed peak 
}
\description{
do the distinction between fishing/non-fishing based on (semi)automatic detection 
of the fishing speed peak to account for the skipper/vessel effect.
This is done per vessel and independently for each gear
of the vessel. A segmented regression is performed on the cumulative speed distribution
to identify the speed breakpoints bounding the fishing peak i.e. the first peak.
The method returns tacsat (a subset if vessels informed) with the SI_STATE informed.
}
\usage{segmentTacsatSpeed(tacsat,vessels,general)}
\arguments{
  \item{tacsat}{tacsat with LE_GEAR informed for each ping}
  \item{general}{a list of general settings e.g. the output path for saving the plot}
  \item{...}{vessels, a vector of vessel identifiers VE_REF that can be found
             in tacsat.
             force.lower.bound might also be used to fix the lower
             speed at some value for gillnet and seine gear types}
}
\details{
A (quick and dirty) vessel and gear-specific method has been suggested by Bastardie et al. 2010 
for dividing the fishing positions from the non-fishing positions (or 'steaming') 
based on vessel speed profile. Data from trips with observers on board recording 
vessel movement with a high spatial resolution (e.g. a position recorded every minute) 
demonstrated that the vessel speed frequency histogram is usually bimodal 
for commercial trawlers, where the first low speed peak likely corresponds to 
fishing activity, while the upper peak corresponds to faster movements, e.g. steaming 
between fishing grounds and harbours. Low speeds can also correspond, for instance, 
to a regulated speed leaving a harbour. 

Speed boundaries were accordingly determined by applying a segmented regression to the 
cumulative distribution of vessel speeds ('calculated' or 'observed', see the argument of the function)
as a practical way to automatically detect 
break points encapsulating the first peak of the speed frequency histogram. 
The computation is preferably done on calculated speed instead of instantaneous (observed) speeds to avoid 
possible non-representative records of speed (discrete records usually with large time-span between). 
Peak detection is done for each vessel and for each of its gear types to account for 
individual skipper behaviour and gear-specific constraints on fishing activity. 
The present methodology is however not well suited for the particular gillnet and seine activities because 
no peak are found from observations of speed histograms. For these fishing activities, 
only the upper boundary of the automatic peak detection should be retained, 
assuming that fishing start at a given value in speed e.g. 0.5 nm.
}
\references{Bastardie et al. 2010. Fisheries Research}
\author{Francois Bastardie}
\seealso{\code{mergeEflalo2Pings}}
\examples{

\dontrun{
data(tacsat) 
# fake here, but an informed LE_GEAR is required
tacsat$LE_GEAR <-"OTB"
tacsat <- segmentTacsatSpeed (
              tacsat=tacsat,
              vessels="35",
              force.lower.bound= 0.5,
              gears.to.force= c('GNS','GND','GNC','GNF','GTR','GTN','GEN','GN','SDN','SSC'),
              general=
                list(
                  output.path=file.path('C:','output'),
                  visual.check=TRUE, speed="calculated"
                  )
              )
}

}