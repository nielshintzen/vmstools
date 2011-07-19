pointOnLand <- function(tacsat,lands){
      if(class(lands) != "SpatialPolygons") stop("'lands' must be specified as class 'SpatialPolygons'")

      rangeLands  <- lapply(lands@polygons[[1]]@Polygons,function(x){c(range(coordinates(x)[,1],na.rm=T),range(coordinates(x)[,2],na.rm=T))})
      lon         <- tacsat$SI_LONG
      lat         <- tacsat$SI_LATI
      res         <- rep(0,length(lon)); totres <- res

      #-Scan over the lands which are given in a list
      for(iLand in 1:length(lands@polygons[[1]]@Polygons)){
        idx     <- which(lon >= rangeLands[[iLand]][1] & lon <= rangeLands[[iLand]][2] & lat >= rangeLands[[iLand]][3] & lat <= rangeLands[[iLand]][4])
        if(length(idx) > 0){
          #- Points on land are 1 or 2 and not is 0
          res[idx]<- point.in.polygon(lon[idx],lat[idx],
                                      coordinates(lands@polygons[[1]]@Polygons[[iLand]])[,1],
                                      coordinates(lands@polygons[[1]]@Polygons[[iLand]])[,2])
          res[idx][which(res[idx] == 2)] <- 0 #no 2 possible
          totres  <- pmax(totres,res)
        }
      }
return(totres)}


