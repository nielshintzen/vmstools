ICESarea <- function(tacsat,areas,proj4string=NULL){
      if(class(areas) != "SpatialPolygons") stop("'areas' must be specified as class 'SpatialPolygons'")

      totres      <- rep(0,length(tacsat$SI_LONG))
      nms         <- unlist(lapply(areas@polygons,function(x){return(x@ID)}))

      if(is.null(proj4string)==T){
        #No projection string used
        if(is.na(proj4string(areas))==F) stop("Projection defined for areas, use proj4string argument in function")
        spPoint           <- SpatialPoints(data.frame(x=tacsat$SI_LONG,y=tacsat$SI_LATI))
        idx               <- over(spPoint,areas)
        totres            <- nms[idx]

      } else {
          #Use projection string
          proj4string(lands)<- proj4string
          spPoint           <- SpatialPoints(data.frame(x=tacsat$SI_LONG,y=tacsat$SI_LATI),proj4string=proj4string)
          idx               <- over(spPoint,areas)
          totres            <- nms[idx]
        }
return(totres)}

