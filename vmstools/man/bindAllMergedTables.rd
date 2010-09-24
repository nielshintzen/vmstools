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
 
 data(eflalo2)
  data(tacsat)
  data(harbours)
  
  # need the tacsat+ format
  tacsat$SI_HARB <- NA
  library(doBy)
  inHarb <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,harbours=harbours,30)
  tacsat$SI_FT <- 1 # init
  idx <- which(inHarb==0)
  tacsat[idx,"SI_FT"] <- cumsum(inHarb) [idx] # add a SI_FT index
  tacsat$SI_STATE <- 2 # init (1: fishing; 2: steaming)
  tacsat$SI_STATE [(tacsat$SI_SP>4 & tacsat$SI_SP<8)] <-1 # fake speed rule for fishing state

  # debug: change funny names of vesselid
  eflalo2$VE_REF <- matrix(unlist(strsplit(as.character(eflalo2$VE_REF),":")),ncol=2,byrow=T)[,2]
  tacsat$VE_REF <- matrix(unlist(strsplit(as.character(tacsat$VE_REF),":")),ncol=2,byrow=T)[,2]

  
  # TEST FOR A GIVEN SET OF VESSELS
  mergeTacsat2EflaloAndDispatchLandingsAtThePingScale (logbooks=eflalo2, tacsat=tacsat, a.vesselid=c("35", "1518"),
                                                             general=list(output.path=file.path("C:","output"),
                                                                            a.year=2009, visual.check=TRUE))
  # ...OR APPLY FOR ALL VESSELS IN eflalo2
  mergeTacsat2EflaloAndDispatchLandingsAtThePingScale (logbooks=eflalo2, tacsat=tacsat,
                                                             general=list(output.path=file.path("C:","output"),
                                                                            a.year=2009, visual.check=TRUE))
  gc(reset=TRUE)

  # load the merged output table for one vessel
  load(file.path("C:","output","merged_35_2009"))
  
  # ...or bind all vessels
  bindAllMergedTables (vessels=c("35", "1518"), species.to.merge=character(), what=character(), 
                      folder = file.path("C:","output"))
 
               
  # map landing of POK
  df1<- merged[,colnames(merged)%in% c("SI_LATI","SI_LONG","LE_KG_POK")]
  df1$SI_LONG <-as.numeric(as.character(df1$SI_LONG))
  df1$SI_LATI <-as.numeric(as.character(df1$SI_LATI))
  vmsGridCreate(df1,nameLon="SI_LONG",nameLat="SI_LATI",cellsizeX =0.05,cellsizeY =0.05)

  # remove steaming points before gridding!
  df2<-df1[-which(df1$LE_KG_POK==0),]
  df3<-df2[-which(is.na(df2$ LE_KG_POK)),]
  vmsGridCreate(df3,nameLon="SI_LONG",nameLat="SI_LATI",cellsizeX =0.05,cellsizeY =0.05)

}
\keyword{datasets}
