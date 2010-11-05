\name{bindAllMergedTables}
\alias{bindAllMergedTables}
\docType{data}
\title{
bind all individual merged tables into a big one 
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(bindAllMergedTables)}
\format{
  The format is:
 chr "bindAllMergedTables"
}
\details{
Possibly, merge irrelevant species into an other category to reduce the size of the output object.
Possibly, keep only the weight, or the value per species. Both by default, but possible out of memory crash. 
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\author{
Francois Bastardie
}
\examples{
 
\dontrun{
  data(eflalo2)
  data(tacsat)
  data(euharbours)
  # add missing harbours? (still to be fix...)
  #euharbours <- c(euharbours, list(a.harbour1=data.frame(lon='10',lat='10')))
  #euharbours <- c(euharbours, list(a.harbour2=data.frame(,lon='1',lat='1')))

 
  library(doBy)
  tacsat$SI_HARB <- NA
  inHarb <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,harbours=euharbours,30)
  tacsat$SI_FT <- 1 # init
  idx <- which(inHarb==0)
  tacsat[idx,"SI_FT"] <- cumsum(inHarb) [idx] # add a SI_FT index
  tacsat <- tacsat[which(inHarb==0),] # keep out of harbour points only
  tacsat$SI_STATE <- 2 # init (1: fishing; 2: steaming)
  tacsat$SI_STATE [(tacsat$SI_SP>4 & tacsat$SI_SP<8)] <-1 # fake speed rule for fishing state


                       
  
  # reduce the size of the eflalo data by merging species (e.g. <1 millions euros)
  eflalo <- mergeEflaloSpecies (eflalo2, threshold=1e6) 
  
  # debug
  eflalo2 <- eflalo2[!eflalo2$VE_REF=="NA" &!is.na(eflalo2$VE_REF),]
  
  # TEST FOR A GIVEN SET OF VESSELS
  mergeTacsat2EflaloAndDispatchLandingsAtThePingScale (logbooks=eflalo2, tacsat=tacsat, a.vesselid=c("35", "1518"),
                                                             general=list(output.path=file.path("C:","output"),
                                                                            a.year=2009, visual.check=TRUE,
                                                                             do.wp3=FALSE, speed="segment"))
  # ...OR APPLY FOR ALL VESSELS IN eflalo2
  mergeTacsat2EflaloAndDispatchLandingsAtThePingScale (logbooks=eflalo2, tacsat=tacsat,
                                                             general=list(output.path=file.path("C:","output"),
                                                                            a.year=2009, visual.check=TRUE,
                                                                             do.wp3=FALSE, speed="segment"))
  gc(reset=TRUE)

  # load the merged output table for one vessel
  load(file.path("C:","output","merged_35_2009.RData"))
  
  # ...or bind all vessels
  tmp <- bindAllMergedTables (vessels=c("35", "1518"), species.to.merge=character(), 
                      folder = file.path("C:","output"), all.in.one.table=TRUE)
 
   # ...and load the merged output table for all vessels
  load(file.path("C:","output","all_merged_2009.RData"))
             
  # map landing of sole from all studied vessels
  df1<- all.merged[,colnames(all.merged)%in% c("SI_LATI","SI_LONG","LE_KG_SOL")]
  df1$SI_LONG <-as.numeric(as.character(df1$SI_LONG))
  df1$SI_LATI <-as.numeric(as.character(df1$SI_LATI))
  vmsGridCreate(df1,nameLon="SI_LONG",nameLat="SI_LATI",cellsizeX =0.05,cellsizeY =0.05)

  # remove steaming points before gridding!
  df2<-df1[-which(is.na(df1$LE_KG_SOL)),]
  vmsGridCreate(df2,nameLon="SI_LONG",nameLat="SI_LATI",cellsizeX =0.05,cellsizeY =0.05)


  # CONVERT TO FISHFRAME FORMAT (might take some time running)
  # (by default, will keep all the species in the output table)
  tmp <- bindAllMergedTables (vessels= unique(tacsat$VE_REF), species.to.keep=character(), 
                      folder = file.path("C:","output"), all.in.one.table=FALSE) 
  
  ff <- mergedTable2Fishframe (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="NLD") )

  }

}
\keyword{datasets}
