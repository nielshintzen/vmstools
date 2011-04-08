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

    load("C:\\merging\\redo3merged2010\\all_merged_weight_2010.RData")
    plotTreeMap (all.merged,  gridcell=c(0.1,0.1), gear="OTB", 
                xlim= c(-1,17), ylim= c(52,62), acolors=rainbow(7),
                  species.to.keep= c("LE_KG_COD","LE_KG_NEP", "LE_KG_PLE", "LE_KG_SOL") ) 

}

}