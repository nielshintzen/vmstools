surface <- function(lon1,lat1,lon2,lat2,res=10){

            x1      <- lon1
            x2      <- lon2
            y1      <- lat1
            y2      <- lat2
            
            if(length(x1) != length(x2) & length(y1) != length(y2) & length(x1) != length(y1)) stop("dimensions of longitudes and latitudes differ")
            
            if(length(x1) == 1){
              sy      <- seq(y1,y2,length.out=res)
              height  <- distance(rep(x1,res-1),sy[2:res],rep(x1,res-1),sy[1:(res-1)])
              base    <- distance(rep(x1,res),sy,rep(x2,res),sy)
              base1   <- base[1:(res-1)]
              base2   <- base[2:res]
              surface <- sum(height * (base1 + base2) / 2)
            }
            if(length(x1) > 1){
              sy      <- mapply(seq,y1,y2,length.out=res)
              sx1     <- mapply(rep,x1,res)
              sx2     <- mapply(rep,x2,res)
              height  <- distance(sx1[1:(res-1),],sy[2:res,],sx1[1:(res-1),],sy[1:(res-1),])
              base    <- distance(sx1,sy,sx2,sy)
              base1   <- base[1:(res-1),]
              base2   <- base[2:res,]
              surface <- apply(height * (base1 + base2) / 2,2,sum)
            }

         return(surface)}