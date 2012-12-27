\name{plotTreeMap}
\alias{plotTreeMap}
\title{
Plot a squarified treemap of landing proportion per cell.
}
\description{
Plot a squarified treemap of landing proportion per cell from the merged VMS/logbooks tables.
}
\usage{
 plotTreeMap (x,  gridcell, gear, 
                xlim, ylim, acolors, species.to.keep) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a merged data.frame such as the one produced by mergeEflalo2Pings()
}
  \item{gridcell}{
  grid cell resolution in degree
}
  \item{gear}{
  subset the data according to a(set of) gear(s) 
}
  \item{xlim}{
  longitude range
}
  \item{ylim}{
  latitude range
}
  \item{acolors}{
  should be at least same of length of the species.to.keep
}
  \item{species.to.keep}{
  set of eflalo names for species (the same for plotting cash value not implemented yet.)
}
}
\details{
to be further custumized...
}
\value{
plot a graph
}
\references{EU Lot 2 project}
\author{Francois Bastardie}
\examples{

 \dontrun{

  
  # do mergeEflalo2Pings()....
  # ...and load the merged output table for all vessels
  load(file.path("C:","output","all_merged__1800.RData"))


  graphics.off()
  df1 <- all.merged[, c("LE_MET_level6", "LE_GEAR","SI_STATE", "SI_LATI",
                        "SI_LONG","LE_KG_COD","LE_KG_PLE")]
  df1$SI_LONG <- anf(df1$SI_LONG)
  df1$SI_LATI <- anf(df1$SI_LATI)
  df1 <-   df1[ !is.na(df1$SI_LATI),]
  df1 <-   df1[ !is.na(df1$SI_LONG),]
  df1$LE_MET  <- df1$LE_MET_level6


  # call to plotTreeMap()
  a.metier <- "TBB_DEF_70-99_0_0"
  a.gear   <-"TBB"
  plotTreeMap (df1[df1$LE_MET_level6 %in% a.metier & df1$SI_STATE==1,],  gridcell=c(0.1,0.05), gear=a.gear,
                xlim= c(3,6), ylim= c(50,54), acolors=rainbow(7),
                  species.to.keep= c("LE_KG_COD","LE_KG_PLE") )
  mtext("Latitude", 2, 3) ;  mtext("Longitude",1, 2)

}
                                 
}