  # Author: F.Bastardie


#' generate some fishing effort maps from the merged output table
#' 
#' generate some effort maps (jpeg files) from the merged table per year, per
#' metier (level6), per metier-quarter all together
#' 
#' using some arguments informing mapGrid(). Performing automatic creation of a
#' hierarchy of output folders.  see also pings2LandingsMaps().
#' 
#' @param all.merged the data.frame from merged VMS with logbooks
#' @param output a path for the output folder to be created
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
#'  years= c('1800')
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
#' }
#' 
#' 
#' @export pings2EffortMaps
  pings2EffortMaps <-
     function (all.merged,  output= file.path("C:","VMSanalysis", "FemernBelt"),
        cellsizeX =0.05, cellsizeY =0.05, we=9.8, ea=12.7, no=55.2, so=54.0,
         breaks0=c(0,25, 50,100,200,400,800,1600, 3200,6400,12800, 100000)){

    # add a quarter
    cat(paste("add a quarter...can be deadly long!", "\n"))
    if(!"quarter" %in% colnames(all.merged))
          all.merged$quarter <- factor(substr(quarters(as.POSIXct(all.merged$SI_DATE)),2,2))

    # create the general folder
    dir.create(file.path(output), recursive = TRUE )
    dir.create(file.path(output, "jpegEffort"))

    # detect the year
    a.year <- format(strptime( paste(all.merged$SI_DATE[1]) , tz='GMT',  "%e/%m/%Y" ), "%Y")


    # in overall
    cat(paste("overall...", "\n"))
    df1 <- all.merged[all.merged$SI_STATE==1, c("SI_LATI","SI_LONG","LE_EFF_VMS")]
    df1$SI_LATI <- anf(df1$SI_LATI ) # debug...
    df1$SI_LONG <- anf(df1$SI_LONG ) # debug...
    df1$LE_EFF_VMS <-anf(df1$LE_EFF_VMS) # debug...
    df1 <-   df1[ !is.na(df1$SI_LATI),]
    df1$LE_EFF_VMS <-   anf(df1$LE_EFF_VMS)/ 60 # in hours
    vmsGridCreate(df1, nameVarToSum="LE_EFF_VMS",  numCats=10,  plotPoints =FALSE, legendtitle="fishing (hours)",
          colLand="darkolivegreen4", addICESgrid=TRUE,
            nameLon="SI_LONG", nameLat="SI_LATI", cellsizeX =cellsizeX, cellsizeY =cellsizeY, we=we, ea=ea, no=no, so=so,
             breaks0=breaks0
           , legendncol=2)  
    dir.create(file.path(output, "jpegEffort", "overall"))
    savePlot(filename = file.path("C:","VMSanalysis", "FemernBelt", "jpegEffort", "overall",
                            paste("map_effort_hours_merged_vessels_",a.year,".jpeg",sep="")),type ="jpeg")
    dev.off()


      # per quarter
    cat(paste("per quarter...", "\n"))
   for (a.quarter in levels(all.merged$quarter) ){
     df1 <- all.merged[all.merged$quarter==a.quarter  &
                all.merged$SI_STATE==1,  c("SI_LATI","SI_LONG","LE_EFF_VMS")]
     df1$SI_LATI <- anf(df1$SI_LATI ) # debug...
     df1$SI_LONG <- anf(df1$SI_LONG ) # debug...
     df1$LE_EFF_VMS <-anf(df1$LE_EFF_VMS) # debug...
     df1 <-   df1[ !is.na(df1$SI_LATI),]
     df1$LE_EFF_VMS <-   anf(df1$LE_EFF_VMS)/ 60 # in hours
      if(nrow(df1)!=0){
        vmsGridCreate(df1, nameVarToSum="LE_EFF_VMS",  numCats=10,  plotPoints =FALSE, legendtitle="fishing (hours)",
          colLand="darkolivegreen4", addICESgrid=TRUE,
            nameLon="SI_LONG", nameLat="SI_LATI", cellsizeX =cellsizeX, cellsizeY =cellsizeY, we=we, ea=ea, no=no, so=so,
             breaks0=breaks0
           , legendncol=2)  
         title(paste(a.year, "-", a.quarter) )
      # create folder and save
      dir.create(file.path(output, "jpegEffort", "quarter"))
      savePlot(filename = file.path(output, "jpegEffort",  "quarter",
        paste("map_effort_hours_merged_vessels_",a.year, a.quarter,".jpeg",sep="")),type ="jpeg")
      dev.off()
      }
     }


    # per metier
    cat(paste("per metier...", "\n"))
    for (met in levels(all.merged$LE_MET_level6) ){
      df1 <- all.merged[all.merged$LE_MET_level6==met &
                all.merged$SI_STATE==1, c("SI_LATI","SI_LONG","LE_EFF_VMS")]
      df1$SI_LATI <- anf(df1$SI_LATI ) # debug...
      df1$SI_LONG <- anf(df1$SI_LONG ) # debug...
      df1$LE_EFF_VMS <-anf(df1$LE_EFF_VMS) # debug...
      df1 <-   df1[ !is.na(df1$SI_LATI),]
      df1$LE_EFF_VMS <-   anf(df1$LE_EFF_VMS)/ 60 # in hours
     if(nrow(df1)!=0){
       vmsGridCreate(df1, nameVarToSum="LE_EFF_VMS",  numCats=10,  plotPoints =FALSE, legendtitle="fishing (hours)",
          colLand="darkolivegreen4", addICESgrid=TRUE,
            nameLon="SI_LONG", nameLat="SI_LATI", cellsizeX =cellsizeX, cellsizeY =cellsizeY, we=we, ea=ea, no=no, so=so,
             breaks0=breaks0
           , legendncol=2)  
          title(paste(met) )
      a.met <- gsub(">=", replacement="o",met) # debug
      a.met <- gsub("<", replacement="",a.met)  # debug
      a.met <- gsub(">", replacement="o",a.met)  # debug
      # create folders and save
      dir.create(file.path(output, "jpegEffort", "overall", a.met))
      savePlot(filename = file.path(output, "jpegEffort","overall", a.met,
                            paste("map_effort_hours_merged_vessels_", a.met,"_",a.year,".jpeg",sep="")),type ="jpeg")
      dev.off()
      }
   }

   # per metier, quarter
    cat(paste("per metier, per quarter...", "\n"))
   for (met in levels(all.merged$LE_MET_level6) ){
     for (a.quarter in levels(all.merged$quarter) ){

      df1 <- all.merged[all.merged$LE_MET_level6==met &
                all.merged$quarter==a.quarter  &
                all.merged$SI_STATE==1,   c("SI_LATI","SI_LONG","LE_EFF_VMS")]
      df1$SI_LATI <- anf(df1$SI_LATI ) # debug...
      df1$SI_LONG <- anf(df1$SI_LONG ) # debug...
      df1$LE_EFF_VMS <-anf(df1$LE_EFF_VMS) # debug...
      df1 <-   df1[ !is.na(df1$SI_LATI),]
      df1$LE_EFF_VMS <-   anf(df1$LE_EFF_VMS)/ 60 # in hours
      if(nrow(df1)!=0){
          vmsGridCreate(df1, nameVarToSum="LE_EFF_VMS",  numCats=10,  plotPoints =FALSE, legendtitle="fishing (hours)",
          colLand="darkolivegreen4", addICESgrid=TRUE,
            nameLon="SI_LONG", nameLat="SI_LATI", cellsizeX =cellsizeX, cellsizeY =cellsizeY, we=we, ea=ea, no=no, so=so,
             breaks0=breaks0
           , legendncol=2)  
      title(paste(met, "-", a.year, "-", a.quarter) )
      a.met <- gsub(">=", replacement="o",met) # debug
      a.met <- gsub("<", replacement="",a.met) # debug
      a.met <- gsub(">", replacement="o",a.met)# debug
      dir.create(file.path(output, "jpegEffort", "quarter", a.met))
      # save
      savePlot(filename = file.path(output, "jpegEffort",  "quarter",  a.met,
        paste("map_effort_hours_merged_vessels_", a.met,"_",a.year, a.quarter,".jpeg",sep="")),type ="jpeg")
      dev.off()
      }
     }}

 return()
 }
 
 
#  pings2EffortMaps (all.merged=all.merged,  output= file.path("C:","VMSanalysis", "FemernBelt"),
#        cellsizeX =0.05, cellsizeY =0.05, we=9.8, ea=12.7, no=55.2, so=54.0,
#         breaks0=c(0,25, 50,100,200,400,800,1600, 3200,6400,12800, 100000))
