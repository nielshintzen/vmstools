  # Author: F.Bastardie


#' generate some landings maps from the merged output table
#' 
#' generate landings maps from the merged table per year, per metier (DCF
#' Level6), per metier-quarter all together
#' 
#' using some arguments informing mapGrid(). Performing automatic creation of a
#' hierarchy of output folders.
#' 
#' @param all.merged the data.frame from merged VMS to logbooks
#' @param sp sp name (FAO code)
#' @param output a path for the output to be located
#' @param cellsizeX additional arguments for mapGrid
#' @param cellsizeY additional arguments for mapGrid
#' @param we additional arguments for mapGrid
#' @param ea additional arguments for mapGrid
#' @param no additional arguments for mapGrid
#' @param so additional arguments for mapGrid
#' @param breaks0 additional arguments for mapGrid
#' @author Francois Bastardie
#' @examples
#' 
#'  
#' \dontrun{
#'  
#'   years= c('1800')
#'   
#'   for(a.year in years){
#'   
#'   graphics.off()
#'   
#'   # ...and load the merged output table for all vessels
#'    load(file.path("C:","output",paste("all_merged__",a.year,".RData",sep='')))
#' 
#'   # generate the effort maps (jpeg files) and store in a hierarchy of folders
#'    pings2EffortMaps (all.merged=all.merged,  output=
#'     file.path("C:","output"),
#'     cellsizeX =0.1, cellsizeY =0.05, we = 3, ea = 6, so = 50, no = 54,
#'      breaks0=c(0,25, 50,100,200,400,800,1600, 3200,6400,12800, 100000))
#' 
#'   # generate the landings maps (jpeg) and store in a hierarchy of folders
#'    pings2LandingsMaps (all.merged=all.merged, sp="LE_KG_COD",
#'          output= file.path("C:","output"),
#'          cellsizeX =0.1, cellsizeY =0.05, we = 3, ea = 6, so = 50, no = 54,
#'          breaks0= c(0,100, 100*(2^1),100*(2^2),100*(2^3),100*(2^4),
#'          100*(2^5),100*(2^6), 100*(2^7),100*(2^8),100*(2^9), 10000000))
#' 
#'   }
#'   
#'  # create some animations from the sequential jpeg files...
#'  # (need the 'animation' R package and ImageMagik installed)
#'  landingsMaps2GIFanim(idir = file.path("C:","output"), spp=c("COD") )
#'     
#'     
#' }
#' 
#' @export pings2LandingsMaps
  pings2LandingsMaps <-
     function (all.merged, sp="LE_KG_COD", output= file.path("C:","VMSanalysis", "FemernBelt"),
        cellsizeX =0.05, cellsizeY =0.05, we=9.8, ea=12.7, no=55.2, so=54.0,
          breaks0= c(0,100, 100*(2^1),100*(2^2),100*(2^3),100*(2^4),100*(2^5),100*(2^6), 100*(2^7),100*(2^8),100*(2^9), 10000000)){

    if(!"quarter" %in% colnames(all.merged))
          all.merged$quarter <- factor(substr(quarters(as.POSIXct(all.merged$SI_DATE)),2,2))

    # create the general folder
    dir.create(file.path(output), recursive = TRUE )
    dir.create(file.path(output, "jpegLandings"))

    # detect the year
    a.year <- format(strptime( paste(all.merged$SI_DATE[1]) , tz='GMT',  "%e/%m/%Y" ), "%Y")

    if (any(grep("EURO", sp)>0)) what <- "value"
    if (any(grep("KG", sp)>0))  what <- "weight"
    if(what=="weight") a.unit <- "(KG)"
    if(what=="value") a.unit <- "(EURO)"

    get.sp <- function (nm) unlist(lapply(strsplit(nm, split="_"), function(x) x[3]))
    a.sp <- get.sp(sp)

    df1 <- all.merged[all.merged$SI_STATE==1, c("SI_LATI","SI_LONG",sp)]
    df1$SI_LATI <- anf(df1$SI_LATI ) # debug...
    df1$SI_LONG <- anf(df1$SI_LONG ) # debug...
    df1[,sp] <- replace(df1[,sp], is.na(df1[,sp]) | df1[,sp]<0, 0)
    vmsGridCreate(df1, nameVarToSum=sp,  numCats=10,  plotPoints =FALSE, legendtitle=paste("landings",what, a.unit,sep=' '),
          colLand="darkolivegreen4",  addICESgrid=TRUE,
            nameLon="SI_LONG", nameLat="SI_LATI", cellsizeX =cellsizeX, cellsizeY =cellsizeY, we=we, ea=ea, no=no, so=so,
             breaks0=breaks0, legendncol=2)
    title(a.sp)
          # create folders and save
    dir.create(file.path(output, "jpegLandings", a.sp))
    dir.create(file.path(output, "jpegLandings", a.sp, "overall"))
    dir.create(file.path(output, "jpegLandings", a.sp, "overall",what))
    savePlot(filename = file.path(output, "jpegLandings", a.sp, "overall", what,
                            paste("map_landings_",what,"_merged_vessels_",a.sp,"_",a.year,".jpeg",sep="")),type ="jpeg")

    dev.off()


      # per quarter
   for (a.quarter in levels(all.merged$quarter) ){
     df1 <- all.merged[all.merged$quarter==a.quarter  &
                all.merged$SI_STATE==1,  c("SI_LATI","SI_LONG",sp)]
      df1$SI_LATI <- anf(df1$SI_LATI) # debug...
      df1$SI_LONG <- anf(df1$SI_LONG) # debug...
      df1[,sp] <- replace(df1[,sp], is.na(df1[,sp]) | df1[,sp]<0, 0)
      if(nrow(df1)!=0){
      vmsGridCreate(df1, nameVarToSum=sp,  numCats=10,  plotPoints =FALSE, legendtitle=paste("landings",what,a.unit,sep=' '),
          colLand="darkolivegreen4",  addICESgrid=TRUE,
            nameLon="SI_LONG", nameLat="SI_LATI", cellsizeX =cellsizeX, cellsizeY =cellsizeY, we=we, ea=ea, no=no, so=so,
             breaks0=breaks0,legendncol=2)
         title(paste( a.sp, "-", a.year, "-", a.quarter) )
      # create folder and save
      dir.create(file.path(output, "jpegLandings", a.sp, "overall", what, "quarter"))
      savePlot(filename = file.path(output, "jpegLandings", a.sp, "overall", what, "quarter",
        paste("map_landings_",what,"_merged_vessels_",a.sp,"_overall_",a.year, a.quarter,".jpeg",sep="")),type ="jpeg")
      dev.off()
      }
     }


    # per metier
    for (met in levels(all.merged$LE_MET_level6) ){
      df1 <- all.merged[all.merged$LE_MET_level6==met &
                all.merged$SI_STATE==1, c("SI_LATI","SI_LONG",sp)]
      df1$SI_LATI <- anf(df1$SI_LATI) # debug...
      df1$SI_LONG <- anf(df1$SI_LONG) # debug...
      df1[,sp] <- replace(df1[,sp], is.na(df1[,sp]) | df1[,sp]<0, 0)
      if(nrow(df1)!=0){
      vmsGridCreate(df1, nameVarToSum=sp,  numCats=10,  plotPoints =FALSE, legendtitle=paste("landings",what,a.unit,sep=' '),
          colLand="darkolivegreen4",  addICESgrid=TRUE,
            nameLon="SI_LONG", nameLat="SI_LATI",cellsizeX =cellsizeX, cellsizeY =cellsizeY, we=we, ea=ea, no=no, so=so,
             breaks0=breaks0,legendncol=2)
         title(paste(met, "-", a.sp) )
      a.met <- gsub(">=", replacement="o",met) # debug
      a.met <- gsub("<", replacement="",a.met)  # debug
      a.met <- gsub(">", replacement="o",a.met)  # debug
      # create folders and save
      dir.create(file.path(output, "jpegLandings", a.sp, a.met))
      dir.create(file.path(output, "jpegLandings", a.sp, a.met, what))
      savePlot(filename = file.path(output, "jpegLandings", a.sp, a.met, what,
                            paste("map_landings_",what,"_merged_vessels_",a.sp,"_", a.met,"_",a.year,".jpeg",sep="")),type ="jpeg")
      dev.off()
      }
   }

   # per metier, quarter
   for (met in levels(all.merged$LE_MET_level6) ){
     for (a.quarter in levels(all.merged$quarter) ){

      df1 <- all.merged[all.merged$LE_MET_level6==met &
                all.merged$quarter==a.quarter  &
                all.merged$SI_STATE==1,  c("SI_LATI","SI_LONG",sp)]
      df1$SI_LATI <- anf(df1$SI_LATI) # debug...
      df1$SI_LONG <- anf(df1$SI_LONG) # debug...
      df1[,sp] <- replace(df1[,sp], is.na(df1[,sp]) | df1[,sp]<0, 0)
      if(nrow(df1)!=0){
      vmsGridCreate(df1, nameVarToSum=sp,  numCats=10,  plotPoints =FALSE, legendtitle=paste("landings",what,a.unit,sep=' '),
          colLand="darkolivegreen4",  addICESgrid=TRUE,
            nameLon="SI_LONG", nameLat="SI_LATI", cellsizeX =cellsizeX, cellsizeY =cellsizeY, we=we, ea=ea, no=no, so=so,
             breaks0=breaks0,legendncol=2)
         title(paste(met, "-", a.sp, "-", a.year, "-", a.quarter) )
      a.met <- gsub(">=", replacement="o",met) # debug
      a.met <- gsub("<", replacement="",a.met) # debug
      a.met <- gsub(">", replacement="o",a.met)# debug
      dir.create(file.path(output, "jpegLandings", a.sp, a.met, what, "quarter"))
      # save
      savePlot(filename = file.path(output, "jpegLandings", a.sp, a.met, what, "quarter",
        paste("map_landings_",what,"_merged_vessels_",a.sp,"_", a.met,"_",a.year, a.quarter,".jpeg",sep="")),type ="jpeg")
      dev.off()
      }
     }}

 return()
 }
  #  in value
 # pings2LandingsMaps (all.merged=all.merged, sp="LE_EURO_COD", output= file.path("C:","VMSanalysis", "FemernBelt"),
 #                   cellsizeX =0.05, cellsizeY =0.05, we=9.8, ea=12.7, no=55.2, so=54.0, # fehmarn Belt area
 #                       breaks0= c(0,100, 100*(2^1),100*(2^2),100*(2^3),100*(2^4),100*(2^5),100*(2^6), 100*(2^7),100*(2^8),100*(2^9), 10000000)
 #                       )
