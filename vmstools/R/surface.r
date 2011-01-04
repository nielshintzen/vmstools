surface <- function(vmsGrid,method="Trapezoid",includeNA=T){  #Methods can be "Trapezoid" and "UTM"
        if(class(vmsGrid) != 'SpatialGridDataFrame') stop("class of vmsGrid should be SpatialGridDataFrame")
        if(!method %in% c("Trapezoid","UTM")) stop("method not available")
        if (class(vmsGrid)=='SpatialGridDataFrame') # not empty...
          {
             if(method == "Trapezoid"){
               res <- max(vmsGrid@grid@cellsize,na.rm=T)/0.1 * 10  #automatic scaling
               if(res < 3) res <- 3
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
          if(method == "UTM"){
            require(PBSmapping)
            griddims <- summary(vmsGrid)$grid
            sizelon  <- griddims[1,2]
            sizelat  <- griddims[2,2]

            ltCentreCell<-coordinates(vmsGrid)

            for (x in 1:(length(ltCentreCell)/2)){
              if (!is.na(vmsGrid@data[x,2]) | includeNA) {        # speed up the calculation by dropping cells with fishing=NA  /!| only work for DCF5 and DCF6!
                minX<-ltCentreCell[x,1]-sizelon/2
                maxX<-ltCentreCell[x,1]+sizelon/2
                minY<-ltCentreCell[x,2]-sizelat/2
                maxY<-ltCentreCell[x,2]+sizelat/2
                ltX<-c(minX,minX,maxX,maxX)
                ltY<-c(minY,maxY,maxY,minY)
                sourcePoly<-cbind(rep(1,4),seq(1,4),ltX,ltY)
                rownames(sourcePoly)<-seq(1,4)
                colnames(sourcePoly)<-c("PID","POS","X","Y")

                polyArea<-calcArea(as.PolySet(sourcePoly, projection="LL"))
                singleCellArea<-polyArea$area
              } else {singleCellArea<-NA}
              vmsGrid@data$cellArea[x]<-singleCellArea
            }
          }
        }
        return(vmsGrid)}
         
         