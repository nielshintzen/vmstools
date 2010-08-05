calculateCI <- function(intLon
                               ,intLat
                               ,vmsIdx1
                               ,vmsIdx2
                               ,VMS.
                               ,grid
                               ,sPDF
                               ,interpolation
                               ,int
                               ,params){

  res1          <- maxRangeCI(intLon,intLat,an(difftime(VMS.$date[vmsIdx2],VMS.$date[vmsIdx1],units="mins")),
                              c(VMS.$speed[vmsIdx1],VMS.$speed[vmsIdx2]))
    #First find the boundaries of the mills ellipse, thereafter, add a 10% extra margin, based on the minimum or
    # maximum value. In the longitude direction, take the minimum value, and find the according latitude to go from km to degrees
  res2          <- range(res1[[1]][,1],na.rm=T); boundx <- c(min(res2) - res1[[2]]*0.2*km2Degree(min(res2),
                                                             res1[[1]][which(min(res2)==res1[[1]][,1]),2],1),
                                                             max(res2) + res1[[2]]*0.2*km2Degree(max(res2),
                                                             res1[[1]][which(max(res2)==res1[[1]][,1]),2],1))
  res3          <- range(res1[[1]][,2],na.rm=T); boundy <- c(min(res3) - res1[[2]]*0.2*1/111.2,
                                                             max(res3) + res1[[2]]*0.2*1/111.2)
  if(res1[[3]] == 1){
    boundx <- c(min(intLon) - res1[[2]]*0.2*km2Degree(min(intLon),
              intLat[which(min(intLon)==intLon)],1),
                max(intLon) + res1[[2]]*0.2*km2Degree(max(intLon),
                intLat[which(max(intLon)==intLon)],1))
    boundy <- c(min(intLat) - res1[[2]]*0.2*1/111.2,
                max(intLat) + res1[[2]]*0.2*1/111.2)
  }

    #Take the upper right and lower left values as setting the boundings of the smaller matrix
  cc2           <- cbind(boundx,boundy)
  idx           <- getGridIndex(cc2,grid,all.inside=F)
    #If grid is too small, then extend grid to fit
  if(any(is.na(idx))){
    grid                  <- createGrid(xrange=cc2[,"boundx"],yrange=cc2[,"boundy"],grid@cellsize[1],grid@cellsize[2])
    spatialGrid           <- SpatialGrid(grid=grid)
    gridded(spatialGrid) = TRUE
    sP                    <- as(spatialGrid,"SpatialPixels")
    sPDF                  <- as(sP,"SpatialPixelsDataFrame")
    sPDF@data             <- data.frame(rep(0,length(sPDF@grid.index)))
    sPDF@data[,2]         <- 0
    colnames(sPDF@data)   <- c("data","tmpdata")
    idx                   <- getGridIndex(cc2,grid,all.inside=T)
  }
    

    #Work out the other elements of the matrix
  row1          <- min(idx)%/%grid@cells.dim[1]+1;            col1          <- min(idx) - (grid@cells.dim[1]*(row1-1))
  row2          <- max(idx)%/%grid@cells.dim[1]+1;            col2          <- max(idx) - (grid@cells.dim[1]*(row2-1))
  pxheigth      <- abs(row2-row1);                            pxwidth       <- abs(col2-col1)
  bbox          <- matrix(NA,nrow=pxheigth+1,ncol=pxwidth+1); bbox[1,]      <- sort(min(idx) - seq(0,pxwidth,1))
  if(pxheigth > 0) for(i in 1:(pxheigth)) bbox[i+1,] <- bbox[1,] + grid@cells.dim[1]*i
  idx           <- c(bbox)

    #Calculate the distan matrix based on the idx
  distan <- matrix(NA,nrow=dim(bbox)[1],ncol=dim(bbox)[2])

  for (x in 2:length(interpolation[[int]][,1])){
    distan <- pmin(distan,
                   matrix(distance(lon=sPDF@coords[idx,1],lat=sPDF@coords[idx,2],lonRef=interpolation[[int]][x,1],latRef=interpolation[[int]][x,2]),
                   nrow=dim(bbox)[1],ncol=dim(bbox)[2]),na.rm=T)
  }
    #Calculate the distance from begin or endpoint
  begindistan <- matrix(distance(lon=sPDF@coords[idx,1],lat=sPDF@coords[idx,2],lonRef=interpolation[[int]][2,1],latRef=interpolation[[int]][2,2]),
                        nrow=dim(bbox)[1],ncol=dim(bbox)[2])
  enddistan   <- matrix(distance(lon=sPDF@coords[idx,1],lat=sPDF@coords[idx,2],lonRef=interpolation[[int]][length(interpolation[[1]][,1]),1],
                                                                               latRef=interpolation[[int]][length(interpolation[[1]][,2]),2]),
                        nrow=dim(bbox)[1],ncol=dim(bbox)[2])
  linepistan  <- pmin(begindistan, enddistan,na.rm=T)
    #Reset very small numbers to 0 to get highest values at begin and end point
  linepistan[which(linepistan < 1e-6)]  <- 0
  distan[which(distan < 1e-6)]          <- 0; zeroDistan <- which(distan==0)

  CI                              <- c(matrix(N1p0(distan*params$distscale,0,linepistan^params$sigline,0),ncol=dim(distan)[2],nrow=dim(distan)[1]))
  if(max(CI,na.rm=T) < 0.1) warning("Prediction max(tmpnew) is very small")
  if(length(zeroDistan)>0)  CI[zeroDistan]  <- pmax(CI[zeroDistan],1,na.rm=T)
  
  return(list(CI,idx,res1,grid,sPDF,ifelse(exists("sP"),sP,0)))}