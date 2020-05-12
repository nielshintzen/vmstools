  # Author: F.Bastardie
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
