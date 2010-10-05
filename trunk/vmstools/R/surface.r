surface <- function(vmsGrid,res=10){
        if (class(vmsGrid)=='SpatialGridDataFrame') # not empty...
          {
             griddims <- summary(vmsGrid)$grid
             bboxdims <- bbox(vmsGrid)
             stlon    <- bboxdims[1,1]
             stlat    <- bboxdims[2,1]
             enlon    <- bboxdims[1,2]
             enlat    <- bboxdims[2,2]
             sizelon  <- griddims[1,2]
             sizelat  <- griddims[2,2]

             lons     <- seq(stlon,enlon,sizelon)
             lats     <- seq(stlat,enlat,sizelat)

             heights  <- distance(lon=0,lat=stlat,lonRef=0,latRef=stlat+sizelat/(res-1))
             seqlats  <- mapply(seq,lats[1:(length(lats)-1)],lats[2:length(lats)],length.out=res)

             base     <- matrix(mapply(distance,lon=0,lat=c(seqlats),lonRef=sizelon,latRef=c(seqlats)),ncol=res,byrow=T)
             if(dim(base)[1] == 1){
                base1     <- base[1:(res-1)]
                base2     <- base[2:res]
                surface   <- rep(sum(heights * (base1 + base2)/2),each=length(seq(stlon,enlon-sizelon,sizelon)))
             } else {
                base1    <- base[,1:(res-1)]
                base2    <- base[,2:res]
                surface <- rep(apply(heights * (base1 + base2) / 2,1,sum),each=length(seq(stlon,enlon-sizelon,sizelon)))
              }

             vmsGrid@data$cellArea <- rev(surface)
          }
        return(vmsGrid)}
         
         