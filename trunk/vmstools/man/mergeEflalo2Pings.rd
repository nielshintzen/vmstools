\name{mergeEflalo2Pings}
\alias{mergeEflalo2Pings}
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
forced to nearest VMS trips to prevent loss of any reported landings.

After this fist merging, the distinction of fishnig and not-fishing states can be performed
if the 'detectFishing' argument is at TRUE. The current default method ('segment') will run an automatic detection 
of the fishing peak given the gear type in use. After having detected the vessel-specific speed lower and upper boundaries,
one state (i.e. in SI_STATE) is assigned to each ping. 

The second step of the process evaluated the consistency between both types of data sources. Sequentially, a quality
flag (i to iv) was issued for each recorded position depending on the various degrees of matching, from
finer to coarser, if both data have in common: (i) the trip identifier, the same area (ICES
rectangle), and the catch date (ii) only the trip identifier and the area, or (iii) only the trip
identifier. (the flag iv occur when not fishing points is remaining while some landings are still to be allocated
in these very few cases, all steaming points are actually converted to fishing points).
The initial datasets (both VMS and logbooks) were split into sub-blocks according to
this flag. Each pair of sub-blocks was then joined separately and all merged sub-blocks were
bound in one dataset afterwards.

In a final step, declared landings (in weight and value) at the ICES rectangle scale from logbooks
were further allocated at a lower spatial and time scale level with regard to the VMS information. 
This was done by allocating a proportion of the landings to each detected fishing position depending on 
the level of matching between VMS and logbook data as described in step 2. In the case of full matching,
landings allocated to a given position were proportional to the number of detected fishing positions lying 
inside each declared logbook area for a given date, assuming that the total landings declared in 
this area at this date are divided equally among all of these positions. In the particular case of 
the use of several gear types within the same rectangle and also for the same catch date, an equal share of 
landings between gear types was also assumed for this day. In case of partial matching, the 
following procedure was used. First, when the match was correct by area but not by catch date, 
the landings were equally dispatched between all the fishing positions of this trip detected in
this particular ICES rectangle. Second, when the match failed on both catch date and area, 
the landings concerned were equally allocated to all the detected fishing positions of the trip 
for which no previous match has been detected. Finally, residual landings (kept in the final output) 
occurred when the match failed both for the catch date and area, and no fishing position remained.

Warning: However, one should take care that the generated maps related to the origin of landings 
do not entirely reflect the spatial distributions 
of the harvested populations but more accurately should be interpreted as 
the conjunction between the stocks and the fishermen effort distribution. 
The distribution of the landings on discrete positions should also be understood with some cautions (
albeit the number of pings may be high enough to be considered representative on a large part of the total exploited areas)
and possible allocation of landings on locations we know that a stock is not present might potentially occur.
}

\usage{
mergeEflalo2Pings(eflalo, tacsat, general=list(output.path=file.path("C:"),
                     visual.check=TRUE, detectFishing=FALSE, speed="segment", 
                      what.speed="calculated", conserve.all=TRUE,
                       ), ...)
}

\arguments{
  \item{eflalo}{data.frame, eflalo format
}
  \item{tacsat}{data.frame, tacsat format
}
  \item{general}{list, general settings to go through all sub-functions e.g. detectFishing
}
  \item{vessels}{(optional) name of the vessel(s) to process e.g. those in tacsat 
}
}
\details{in tacsat SI_STAT, fishing should be coded 1 and non-fishing should be coded 2.
The method is very sensitive to missing harbours. Some warning messages are then provided to help the user
to detect some repeated positions that are likely to be missing harbours. 
In this last case you should revise the harbour list before processing again.
}

\value{
Nothing is returned but a merged data.frame per vessel in the output folder. 
These data.frame could be later bound into a big one using bindAllMergedTable()
}

\references{Bastardie et al. 2010. Fisheries Research}
\author{Francois Bastardie}
\note{}


\seealso{\code{\link{vmsGridCreate}}}

\examples{

  \dontrun{
  data(eflalo)
  data(tacsat)
  data(euharbours)

  # format
  eflalo <- formatEflalo(eflalo)
  tacsat <- formatTacsat(tacsat)

  # order tacsat chronologically with library(doBy)
  tacsat <- sortTacsat(tacsat)

  # test each ping if in harbour or not
  tacsat$SI_HARB <- NA
  euharbours$Description <- euharbours$harbour
  tacsat$SI_HARB <- pointInHarbour(lon=anf(tacsat$SI_LONG),
                                   lat=anf(tacsat$SI_LATI),
                                   harbours=euharbours,
                                   rowSize=30, returnNames=TRUE)
  inHarb <- tacsat$SI_HARB
  inHarb <- replace(inHarb, !is.na(inHarb), 1)
  inHarb <- replace(inHarb, is.na(inHarb), 0)
  inHarb <- as.numeric(inHarb)

  # assign a trip identifier
  tacsat$SI_FT <- 1 # init
  idx <- which(inHarb==0)
  tacsat[idx,"SI_FT"] <- cumsum(inHarb) [idx] # add a SI_FT index

  # keep 'out of harbour' points only
  # (but keep the departure point and the arrival point lying in the harbour)
  startTrip <- c(diff(tacsat[,"SI_FT"]), 0)
  endTrip   <- c(0, diff(tacsat[,"SI_FT"]))
  tacsat[which(startTrip>0),"SI_FT"]  <-  tacsat[which(startTrip>0)+1,"SI_FT"] 
  tacsat[which(endTrip<0),"SI_FT"]    <-  tacsat[which(endTrip<0)-1,"SI_FT"] 
  tacsat <- tacsat[which(inHarb==0 |  startTrip>0 |  endTrip<0),]


  # assign a state to each ping (here, useless if detectFishing at TRUE)
  tacsat$SI_STATE <- 2 # init (1: fishing; 2: steaming)
  # fake speed rule for fishing state
  tacsat$SI_STATE [(tacsat$SI_SP>4 & tacsat$SI_SP<8)] <-1


  # reduce the size of the eflalo data by merging species
  # (assuming that the other species is coded MZZ), threshold in euros.
  eflalo2 <- poolEflaloSpecies (eflalo, threshold=1e6, code="MZZ")

  # debug if eflalo has not been cleaned earlier
  eflalo <- eflalo[!eflalo$VE_REF=="NA" &!is.na(eflalo$VE_REF),]
  
  # an informed VE_FLT is also required
  if(all(is.na(eflalo$VE_FLT))) eflalo$VE_FLT <- "fleet1"
  
  # possible mis-naming mistakes
    if(!match('LE_MET_level6',colnames(eflalo))>0){
      eflalo$LE_MET_level6 <- eflalo$LE_MET
    }

  # debug
  eflalo <- eflalo[eflalo$LE_MET!="No_logbook6",]


  # TEST FOR A GIVEN SET OF VESSELS
  # (if detect.fishing is true then do also detection of fishing activity
  # e.g. if speed='segment' the segmentTacsatSpeed() automatic detection of fishing states
  # that will overwrite the existing SI_STATE)
  mergeEflalo2Pings (eflalo=eflalo, tacsat=tacsat, vessels=c("738", "804"),
                     general=list(output.path=file.path("C:","output"),
                     visual.check=TRUE, detectFishing=TRUE, speed="segment",
                     what.speed="calculated"))
  # ...OR APPLY FOR ALL VESSELS IN eflalo
  mergeEflalo2Pings (eflalo=eflalo, tacsat=tacsat,
                     general=list(output.path=file.path("C:","output"),
                     visual.check=TRUE, detectFishing=TRUE, speed="segment",
                     what.speed="calculated"))
  gc(reset=TRUE)

  # load the merged output table for one vessel
  load(file.path("C:","output","merged_804_1800.RData"))

  # check the conservation of landings
  sum(tapply(anf(merged$LE_KG_PLE), merged$flag, sum, na.rm=TRUE))
  sum(eflalo[eflalo$VE_REF=="804","LE_KG_PLE"], na.rm=TRUE)


   # ...or bind all vessels (keeping only some given species here)
  bindAllMergedTables (vessels=c("738", "804"), a.year = "1800",
                      species.to.keep=c("PLE","COD"),
                      folder = file.path("C:","output"),
                      all.in.one.table=TRUE)

    # ...and load the merged output table for all vessels
  load(file.path("C:","output","all_merged__1800.RData"))

  # map landing of cod from all studied vessels
  # ( with debugging if tacsat has not been cleaned earlier)
  graphics.off()
  df1<- all.merged[, c("SI_LATI","SI_LONG","LE_KG_COD")]
  df1$SI_LONG <- anf(df1$SI_LONG)
  df1$SI_LATI <- anf(df1$SI_LATI)
  df1 <-   df1[ !is.na(df1$SI_LATI),]
  df1 <-   df1[ !is.na(df1$SI_LONG),]
  vmsGridCreate(df1,nameLon="SI_LONG", nameLat="SI_LATI",
                nameVarToSum = "LE_KG_COD", cellsizeX =0.1,
                cellsizeY =0.05,  legendtitle = "COD landings (kg)")

 # but you need to remove steaming points before gridding!
  df2<-df1[-which(is.na(df1$LE_KG_COD)),]
  vmsGridCreate(df2,nameLon="SI_LONG",nameLat="SI_LATI", we = 3, ea = 6, so = 50, no = 54,
                nameVarToSum = "LE_KG_COD",cellsizeX =0.1,
                cellsizeY =0.05,  legendtitle = "COD landings (kg)", plotPoints =TRUE, 
                breaks0=c(1,2,4,8,16,32,64,100000))



  # CONVERT TO FISHFRAME FORMAT (might take some time running)
  # (by default, this will keep all the species in the output table)
  tmp <- bindAllMergedTables (vessels= unique(tacsat$VE_REF),
                              species.to.keep=character(),
                              folder = file.path("C:","output"),
                              all.in.one.table=FALSE)

  ff  <- pings2Fishframe (general=list(output.path=file.path("C:","output"),
                          a.year=1800, a.country="NLD", degree=0.05 ) )

  }

}

