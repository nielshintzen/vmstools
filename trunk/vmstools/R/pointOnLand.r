pointOnLand <- function(tacsat,lands){

      rangeLands <- lapply(lands,function(x){c(range(x$x,na.rm=T),range(x$y,na.rm=T))})

      lon       <- tacsat$SI_LONG
      lat       <- tacsat$SI_LATI
      res       <- rep(0,length(lon)); totres <- res

      #-Scan over the lands which are given in a list
      for(iLand in 1:length(lands)){
        idx     <- which(lon >= rangeLands[[iLand]][1] & lon <= rangeLands[[iLand]][2] & lat >= rangeLands[[iLand]][3] & lat <= rangeLands[[iLand]][4])
        if(length(idx) > 0){
          #- Points on land are 1 or 2 and not is 0
          res[idx]<- point.in.polygon(lon[idx],lat[idx],lands[[iLand]]$x,lands[[iLand]]$y)
          res[idx][which(res[idx] == 2)] <- 0 #no 2 possible
          totres  <- pmax(totres,res)
        }
      }
return(totres)}


