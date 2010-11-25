\name{mergeTacsat2EflaloAndDispatchLandingsAtThePingScale}
\alias{mergeTacsat2EflaloAndDispatchLandingsAtThePingScale}
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
if the do.wp3 argument is at TRUE. The current default method ('segment') will run an automatic detection 
of the fishing peak given the gear type in use. After having detected the vessel-specific speed lower and upper boundaries,
one state (i.e. in SI_STATE) is assigned to each ping. 

The second step of the process evaluated the consistency between both types of data sources. Sequentially, a quality
flag was issued for each recorded position depending on the various degrees of matching, from
finer to coarser, if both data have in common: (i) the trip identifier, the same area (ICES
rectangle), and the catch date (ii) only the trip identifier and the area, or (iii) only the trip
identifier. The initial datasets (both VMS and logbooks) were split into sub-blocks according to
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


}

\usage{
mergeTacsat2EflaloAndDispatchLandingsAtThePingScale(logbooks, vms, general = general, ...)
}

\arguments{
  \item{logbooks}{data.frame, eflalo format
}
  \item{tacsat}{data.frame, tacsat format
}
  \item{general}{list, general settings to go through all functions
}
  \item{\dots}{character, a.vesselid = name of the vessel(s) in tacsat (and hopefully, in eflalo)
}
}
\details{in tacsat SI_STAT, fishing should be coded 1 and non-fishing should be coded 2.
The method is very sensitive to missing harbours. Some warning messages are then providen to help the user
to detect some repeated positions that are likely to be missing harbours.
}

\value{
Nothing is returned but a merged data.frame per vessel in the output folder}

\references{Bastardie et al. 2010. Fisheries Research}
\author{Francois Bastardie}
\note{}


\seealso{vmsGridCreate}

\examples{

 \dontrun{
  data(eflalo2)
  data(tacsat)
  data(euharbours)
  # add some missing harbours to the list? 
  #euharbours <- c(euharbours, list(a.harbour1=data.frame(lon='10',lat='10', range='3')))
  #euharbours <- c(euharbours, list(a.harbour2=data.frame(,lon='1',lat='1', range='3')))

  
  # create one column 'date.in.R' from pasting date and time 
  # (to order tacsat if needed)
  #ctime <- strptime(  paste(tacsat$SI_DATE, tacsat$SI_TIME) , tz='GMT',   "%e/%m/%Y %H:%M" )
  tacsat <- cbind.data.frame(tacsat, date.in.R=ctime)
     
  # order tacsat
  library(doBy)
  tacsat <- sortTacsat(tacsat)

  # test each ping if in harbour or not
  tacsat$SI_HARB <- NA
  inHarb <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,harbours=euharbours, 30)
  
  # assign a trip identifier
  tacsat$SI_FT <- 1 # init
  idx <- which(inHarb==0)
  tacsat[idx,"SI_FT"] <- cumsum(inHarb) [idx] # add a SI_FT index
  
  # keep out of harbour points only
  tacsat <- tacsat[which(inHarb==0),] 
  
  # assign a state to each ping (start guesses only)
  tacsat$SI_STATE <- 2 # init (1: fishing; 2: steaming)
  tacsat$SI_STATE [(tacsat$SI_SP>4 & tacsat$SI_SP<8)] <-1 # fake speed rule for fishing state
                     
  
  # reduce the size of the eflalo data by merging species (e.g. <1 millions euros)
  # (assuming that the other species is coded MZZ)
  eflalo <- mergeEflaloSpecies (eflalo2, threshold=1e6, code="MZZ") 
  
  # debug
  eflalo2 <- eflalo2[!eflalo2$VE_REF=="NA" &!is.na(eflalo2$VE_REF),]
  
  # TEST FOR A GIVEN SET OF VESSELS
  # (if do.wp3 is at true then do also the automatic detection of fishing states
  # that will overwrite the existing SI_STATE)
  mergeTacsat2EflaloAndDispatchLandingsAtThePingScale (logbooks=eflalo2, tacsat=tacsat, a.vesselid=c("35", "1518"),
                                                             general=list(output.path=file.path("C:","output"),
                                                                             visual.check=TRUE,
                                                                             do.wp3=TRUE, speed="segment"))
  # ...OR APPLY FOR ALL VESSELS IN eflalo2
  mergeTacsat2EflaloAndDispatchLandingsAtThePingScale (logbooks=eflalo2, tacsat=tacsat,
                                                             general=list(output.path=file.path("C:","output"),
                                                                             visual.check=TRUE,
                                                                             do.wp3=FALSE, speed="segment"))
  gc(reset=TRUE)

  # load the merged output table for one vessel
  load(file.path("C:","output","merged_35_2009.RData"))
  
  # ...or bind all vessels (keeping only some given species here)
  bindAllMergedTables (vessels=c("35", "1518"), species.to.keep=c("PLE","COD"), 
                      folder = file.path("C:","output"), all.in.one.table=TRUE)
 
   # ...and load the merged output table for all vessels
  load(file.path("C:","output","all_merged__2009.RData"))
             
  # map landing of cod from all studied vessels
  df1<- all.merged[, c("SI_LATI","SI_LONG","LE_KG_COD")]
  df1$SI_LONG <-as.numeric(as.character(df1$SI_LONG))
  df1$SI_LATI <-as.numeric(as.character(df1$SI_LATI))
  df1 <-   df1[ !is.na(df1$SI_LATI),]
  df1 <-   df1[ !is.na(df1$SI_LONG),]
  vmsGridCreate(df1,nameLon="SI_LONG",nameLat="SI_LATI", nameVarToSum = "LE_KG_COD",
                                 cellsizeX =0.05,cellsizeY =0.05,  legendtitle = "landings (kg)")

  # remove steaming points before gridding!
  df2<-df1[-which(is.na(df1$LE_KG_COD)),]
  vmsGridCreate(df2,nameLon="SI_LONG",nameLat="SI_LATI", nameVarToSum = "LE_KG_COD",
                                cellsizeX =0.05,cellsizeY =0.05,  legendtitle = "landings (kg)", 
                                 breaks0=c(1,2,4,8,16,32,64,100000))



  # CONVERT TO FISHFRAME FORMAT (might take some time running)
  # (by default, this will keep all the species in the output table)
  tmp <- bindAllMergedTables (vessels= unique(tacsat$VE_REF), species.to.keep=character(), 
                      folder = file.path("C:","output"), all.in.one.table=FALSE)
  
  ff <- mergedTable2Fishframe (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="NLD") )

 
 
  # TO DO....
  # Use the interpolation routine to improve the location of the effort
  #interpolations      <- interpolateTacsat( all.merged 
  #                            ,interval=60             
  #                            ,margin=12               
  #                            ,res=100                
  #                            ,method="cHs"           
  #                            ,params=list(fm=0.5,distscale=20,sigline=0.2,st=c(2,6)) 
  #                            ,headingAdjustment=0
  #                            )
  #interpolationsED <- equalDistance(interpolations,res=10)
  # make sure that the 'res' statement in the interpolateTacsat is significantly bigger 
  # than the 'res' statement in the equalDistance function.
   
  # then map again...
  #vmsGridCreate(interpolationsED,nameLon="SI_LONG",nameLat="SI_LATI", 
  #          cellsizeX =0.05, cellsizeY =0.05, legendtitle = "landings (kg)")

  
 
 
  }

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
