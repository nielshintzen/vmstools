 # function to assign a ICES area code to each ping
 # from the ICES_areas shape file
 ICESarea <- function(long, lat,...){
   lstargs <- list(...)
   library(shapefiles)
   library(sp)
   if(length(lstargs$shape.file)==0) { data(ICESareas)
   } else { ICESareas <- read.shapefile(lstargs$shape.file) }
   an <- function(x) as.numeric(as.character(x))
   ICES_area <- rep(NA, length(long))  # init
   for (i in 1:length(ICESareas$shp$shp)){
     print(as.character(ICESareas$dbf$dbf$ICESAREA[i]))
     a.polygon <- cbind(ICESareas$shp$shp[[i]]$points$X, ICESareas$shp$shp[[i]]$points$Y)
     idx.inout <- point.in.polygon(an(long), an(lat), a.polygon[,1], a.polygon[,2])
     ICES_area[idx.inout==1] <- as.character(ICESareas$dbf$dbf$ICESAREA[i])
     }
 return(ICES_area)
 }
 #all.merged$ICES_area <- ICESarea(long= all.merged$SI_LONG, lat= all.merged$SI_LATI,
 #              shape.file=file.path("H:","DIFRES","VMSanalysis","background_map_shape_files","ICES_areas"))
