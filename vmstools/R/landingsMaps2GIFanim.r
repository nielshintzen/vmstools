 # author F.Bastardie
 # automatic creation of .gif animation using library(animation)
 # need to install the ImageMagik freeware


#' generate some gif animations from the landings maps
#' 
#' generate some gif animations from the landings maps after having generated
#' the landings maps from the merged table per year, per metier (level6), per
#' metier-quarter all together
#' 
#' Using some facilities providen by the 'animation' R package. Need ImageMagik
#' installed.
#' 
#' @param idir the full path of the jpegLandings folder
#' @param spp vector of FAO codes for species
#' @author Francois Bastardie
#' @examples
#' 
#'  
#' \dontrun{
#' 
#' 
#' for(a.year in as.character(2005:2009)){
#' 
#'   # ...and load the merged output table for all vessels
#'   load(file.path("C:","output",paste("all_merged__",a.year,".RData",sep='')))
#' 
#'   # generate the maps
#'   pings2LandingsMaps (all.merged=all.merged, sp="LE_EURO_COD",
#'       output= file.path("C:","VMSanalysis", "FemernBelt"),
#'       cellsizeX =0.05, cellsizeY =0.05, we=9.8, ea=12.7, no=55.2, so=54.0,
#'       breaks0= c(0,100, 100*(2^1),100*(2^2),100*(2^3),100*(2^4),
#'       100*(2^5),100*(2^6), 100*(2^7),100*(2^8),100*(2^9), 10000000))
#'   }
#' 
#' # create some animations...
#' # (need the 'animation' R package and ImageMagik installed)
#' landingsMaps2GIFanim(idir = file.path("C:","VMSanalysis","FemernBelt",
#'                      "jpegLandings"), spp=c("COD") )
#' 
#' 
#'   }
#' 
#' 
#' @export landingsMaps2GIFanim
landingsMaps2GIFanim <- function(idir = file.path("C:","VMSanalysis","FemernBelt","jpegLandings"),
                           spp=c("COD","SPR","HER","PLE","FLE","DAB","WHB")){
require(animation)
for(sp in spp){
    for(met in list.files(file.path(idir, sp))){
     for(what in c("weight","value")){
           setwd(file.path(idir, sp, met,what,"quarter"))
           filename <- paste("map_landings_",what,"_merged_vessels_",sp,"_",met,"_",sep='')
           cat(paste(filename,"\n"))
           wildcard = paste(filename, "*.", "jpeg", sep = "")
           im.convert(wildcard, output = file.path(getwd(), paste("anim",".gif", sep='')),  # library(animation)
                         convert = "convert",
                          cmd.fun = system, clean = FALSE)   
    } # end what
  } # end met
 } # end sp
return()
}

 # pings2LandingMaps (all.merged=all.merged, sp="LE_EURO_COD",  output= file.path("C:","VMSanalysis", "FemernBelt"),
 #                   cellsizeX =0.05, cellsizeY =0.05, we=9.8, ea=12.7, no=55.2, so=54.0, # fehmarn Belt area
 #                       breaks0= c(0,100, 100*(2^1),100*(2^2),100*(2^3),100*(2^4),100*(2^5),100*(2^6), 100*(2^7),100*(2^8),100*(2^9), 10000000)
 #                       )
 # pings2LandingMaps (all.merged=all.merged, sp="LE_KG_COD",  output= file.path("C:","VMSanalysis", "FemernBelt"),
 #                   cellsizeX =0.05, cellsizeY =0.05, we=9.8, ea=12.7, no=55.2, so=54.0, # fehmarn Belt area
 #                       breaks0= c(0,100, 100*(2^1),100*(2^2),100*(2^3),100*(2^4),100*(2^5),100*(2^6), 100*(2^7),100*(2^8),100*(2^9), 10000000)
 #                       )
 # landingsMaps2GIFanim(idir = file.path("C:","VMSanalysis","FemernBelt","jpegLandings"), spp=c("COD","SPR") )
