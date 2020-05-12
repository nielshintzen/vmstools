pointOnLand <- function(tacsat,lands,proj4string=NULL){
      if(class(lands) != "SpatialPolygons") stop("'lands' must be specified as class 'SpatialPolygons'")

      totres      <- rep(0,length(tacsat$SI_LONG))

      if(is.null(proj4string)==TRUE){
        #No projection string used
        if(is.na(proj4string(lands))==FALSE) stop("Projection defined for lands, use proj4string argument in function")
        spPoint           <- SpatialPoints(data.frame(x=tacsat$SI_LONG,y=tacsat$SI_LATI))
        idx               <- over(spPoint,lands)
        totres[which(is.na(idx)==FALSE)] <- 1
        totres[which(is.na(idx)==TRUE)] <- 0

      } else {
          #Use projection string
          proj4string(lands)<- proj4string
          spPoint           <- SpatialPoints(data.frame(x=tacsat$SI_LONG,y=tacsat$SI_LATI),proj4string=proj4string)
          idx               <- over(spPoint,lands)
          totres[which(is.na(idx)==FALSE)] <- 1
          totres[which(is.na(idx)==TRUE)] <- 0
        }
return(totres)}


