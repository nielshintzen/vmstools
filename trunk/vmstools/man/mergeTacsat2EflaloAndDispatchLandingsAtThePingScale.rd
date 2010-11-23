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
\details{in tacsat SI_STAT, fishing should be coded 1 and non-fishing should be coded 2}

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

   # order
  library(doBy)
  tacsat <- orderBy (~VE_REF+SI_DATE+SI_TIME, data=tacsat)

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
                                                                            a.year=2009, visual.check=TRUE,
                                                                             do.wp3=TRUE, speed="segment"))
  # ...OR APPLY FOR ALL VESSELS IN eflalo2
  mergeTacsat2EflaloAndDispatchLandingsAtThePingScale (logbooks=eflalo2, tacsat=tacsat,
                                                             general=list(output.path=file.path("C:","output"),
                                                                            a.year=2009, visual.check=TRUE,
                                                                             do.wp3=FALSE, speed="segment"))
  gc(reset=TRUE)

  # load the merged output table for one vessel
  load(file.path("C:","output","merged_35_2009.RData"))
  
  # ...or bind all vessels (keeping only some given species here)
  tmp <- bindAllMergedTables (vessels=c("35", "1518"), species.to.keep=c("PLE","COD"), 
                      folder = file.path("C:","output"), all.in.one.table=TRUE)
 
   # ...and load the merged output table for all vessels
  load(file.path("C:","output","all_merged__2009.RData"))
             
  # map landing of cod from all studied vessels
  df1<- all.merged[,colnames(all.merged)%in% c("SI_LATI","SI_LONG","LE_KG_COD")]
  df1$SI_LONG <-as.numeric(as.character(df1$SI_LONG))
  df1$SI_LATI <-as.numeric(as.character(df1$SI_LATI))
  vmsGridCreate(df1,nameLon="SI_LONG",nameLat="SI_LATI", 
            cellsizeX =0.05, cellsizeY =0.05, legendtitle = "landings (kg)")

  # but need to remove steaming points before gridding!
  df2<-df1[-which(is.na(df1$LE_KG_COD)),]
  vmsGridCreate(df2,nameLon="SI_LONG",nameLat="SI_LATI", 
             cellsizeX =0.05,cellsizeY =0.05, legendtitle = "landings (kg)")


  # CONVERT TO FISHFRAME FORMAT (might take some time running)
  # (by default, this will keep all the species in the output table)
  tmp <- bindAllMergedTables (vessels= unique(tacsat$VE_REF), species.to.keep=character(), 
                      folder = file.path("C:","output"), all.in.one.table=FALSE)
  
  ff <- mergedTable2Fishframe (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="NLD") )

  }

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
