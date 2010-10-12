 TO FISHFRAME
# create the VE table to upload in fishframe
# required: data.table package
# required: the ICES_areas shape file
 mergedTable2FishframeVE <- function (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="DNK",
                              shape.file=file.path("C:","ICES_areas") )){

    for(what in c('value','weight')){
       load(file.path(general$output.path, paste("all_merged_",what,"_",general$a.year,".RData",sep='')))
       nm <- colnames(all.merged)
       if (what=='value') idx.col  <- grep('EURO', nm) # index columns with species
       if (what=='weight') idx.col  <- grep('KG', nm) # index columns with species
       assign(what, apply(all.merged[,idx.col], 1, sum, na.rm=TRUE))
       }
    all.merged$totvalue <- value
    all.merged$totweight <- weight

    an <- function(x) as.numeric(as.character(x))
    all.merged <- all.merged[, !colnames(all.merged) %in% c('fuelcons', 'flag')]
    all.merged <- all.merged[all.merged$SI_STATE==1,] # keep only fishing pings
    nm <- colnames(all.merged)
    all.merged$ICES_area <- ICESarea(long= all.merged$SI_LONG, lat= all.merged$SI_LATI,
                                         shape.file=general$shape.file) # DEADLY SLOW!
    all.merged$c_square  <-CSquare(an(all.merged$SI_LONG), an(all.merged$SI_LATI), degrees=0.05)
    all.merged$month <- factor(format(as.POSIXct(all.merged$SI_DATE), "%m"))  # add month
    all.merged$LE_VMS_EFF <- an(all.merged$LE_EFF_VMS) / 24 # convert in hours
    all.merged <- all.merged[,c("LE_EFF_VMS","KW_HOURS","totvalue", "totweight", "LE_MET_level6","ICES_area","c_square","month")]
    all.merged$c_square <- factor(all.merged$c_square)
    all.merged$ICES_area <- factor(all.merged$ICES_area)

 
    # base::aggregate() replaced by fast grouping using the data.table library
    library(data.table)
    DT <- data.table(all.merged)
    qu = quote(list(sum(an(LE_EFF_VMS)),sum(an(KW_HOURS)),sum(an(totvalue)),sum(an(totweight))))
    ff.ve <- DT[,eval(qu), by=list(c_square,ICES_area, month,LE_MET_level6)]
    colnames(ff.ve) <- c('c_square','ICES_area', 'month','LE_MET_level6','hours','kw_hours', 'totvalue','totweight')

    # additional
    ff.ve$year       <- general$a.year
    ff.ve$country    <- general$a.country
    ff.ve$nationalFAC    <- NA
    ff.ve$recordtype <- "VE"
    ff.ve$quarter <- ff.ve$month  # init
    levels(ff.ve$quarter) <- c(1,1,1,2,2,2,3,3,3,4,4,4)


  #order colums
  ff.ve <- as.data.frame(ff.ve)[, c('recordtype','country','year','quarter', 'month', 
               'ICES_area','c_square', 'nationalFAC', 'LE_MET_level6',
                 'hours','kw_hours','totweight','totvalue')]
  # save
  write.table(ff.ve, file=file.path(general$output.path, 
         paste("ff_ve_", general$a.year, ".txt", sep='')),  dec=".", sep=";", quote=FALSE, row.names=FALSE)

  return(ff.ve)
 }

 # function to assign a ICES area code to each ping 
 # from the ICES_areas shape file
 ICESarea <- function(long, lat,
       shape.file=file.path("C:","ICES_areas")){
   library(shapefiles) 
   library(sp)
   an <- function(x) as.numeric(as.character(x))
   ICES_area <- rep(NA, length(long))  # init
   a.shape <- read.shapefile(shape.file)  
   for (i in 1:length(a.shape$shp$shp)){       
     print(as.character(a.shape$dbf$dbf$ICESAREA[i]))
     a.polygon <- cbind(a.shape$shp$shp[[i]]$points$X,a.shape$shp$shp[[i]]$points$Y)  
     idx.inout <- point.in.polygon(an(long), an(lat), a.polygon[,1], a.polygon[,2]) 
     ICES_area[idx.inout==1] <- as.character(a.shape$dbf$dbf$ICESAREA[i])
     }           
 return(ICES_area)
 }
 
  # call
 #mergedTable2FishframeVE (general=list(output.path=file.path("H:","DIFRES","VMSanalysis","results_merged","DKWaters"),
 #     a.year=2009, a.country="DNK",shape.file=file.path("H:","DIFRES","VMSanalysis","background_map_shape_files","ICES_areas")) )
