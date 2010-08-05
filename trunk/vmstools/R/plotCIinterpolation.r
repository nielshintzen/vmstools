`plotCIinterpolation` <-
function(interpolation
                                             ,VMS=NULL
                                             ,grid
                                             ,plot=F
                                             ,params=list(0.2,20,0.2)  #Specify the three parameters: fm, distscale, sigline.
                                             ,plotPoint=F
                                             ){

VMS. <- VMS
                                             
  #Define the single points, and the interpolated points
pPoints       <- seq(1,dim(VMS.)[1],1)
pInterOnly    <- unique(unlist(lapply(interpolation,function(x){return(x[1,])})))
pPointsOnly   <- which(is.na(pmatch(pPoints,pInterOnly))==T)

  #Create SpatialGridDataFrame, the object where we can store the impact data
spatialGrid   <- SpatialGrid(grid=grid)
gridded(spatialGrid) = TRUE
sP            <- as(spatialGrid,"SpatialPixels")
sPDF          <- as(sP,"SpatialPixelsDataFrame")
sPDF@data     <- data.frame(rep(0,length(sPDF@grid.index)))
sPDF@data[,2] <- 0
colnames(sPDF@data) <- c("data","tmpdata")

  #Create object to count the difftime to give an estimate for the points-only difftime
timeDiff      <- numeric()

  #Create a matrix with all the inter-only VMS datapoints
cc            <- na.omit(as.matrix(cbind(VMS.$declon[pInterOnly],VMS.$declat[pInterOnly])))

  #Create smaller matrix to work from for a single interpolation
for(int in 1:length(interpolation)){
  iPV1   <- interpolation[[int]][1,1]; iPV2   <- interpolation[[int]][1,2]
  iP1    <- which(iPV1 == pInterOnly); iP2    <- which(iPV2 == pInterOnly)

    #Count the difftime to give an estimate for the points-only difftime
  timeDiff[int] <- an(difftime(VMS.$date[iPV2],VMS.$date[iPV1],units="mins"))
    #Calculate the confidence interval
  resCI <- calculateCI(cc[c(iP1,iP2),1],cc[c(iP1,iP2),2],iPV1,iPV2,VMS.,grid,sPDF,interpolation,int,params)
  CI    <- resCI[[1]]
  idx   <- resCI[[2]]
  grid  <- resCI[[4]]
  sPDF  <- resCI[[5]]
  if(class(resCI[[6]])!="numeric") sP <- resCI[[6]] 

    #Reset the new data set at 0
  sPDF@data[,2]                   <- data.frame(rep(0,length(sPDF@grid.index)))
  if(int == 1){ sPDF@data[idx,1]  <- CI; idxmin1 <- idx; sPDF@data[idx,3] <- CI
  } else {      sPDF@data[idx,2]  <- CI}

    #Determine if the interpolation is overlapping with the preceeding point, if so, compensate for that
  if(int != 1){
      #Three data storages in sPDF@data
        #index 1: final data with all interpolated tracks
        #index 2: temporal data with last interpolated track
        #index 3: temporal data with last - 1 interpolated track

    if(interpolation[[int]][1,1] == (interpolation[[int-1]][1,2])){
      matching                <- idx[which(idx%in%idxmin1==T)]
      nomatching              <- idx[which(idx%in%idxmin1==F)]
      if(length(matching) > 0)    sPDF@data[matching,1]   <- pmax(sPDF@data[matching,3],sPDF@data[matching,2],na.rm=T) +
                                                             sPDF@data[matching,1] - sPDF@data[matching,3]
      if(length(nomatching) > 0)  sPDF@data[nomatching,1] <- sPDF@data[nomatching,2] + sPDF@data[nomatching,1]
      sPDF@data[idx,3]        <- sPDF@data[idx,2]
      idxmin1                 <- idx
    } else {
        sPDF@data[idx,1]      <- CI + sPDF@data[idx,1]
        sPDF@data[idx,3]      <- sPDF@data[idx,2]
        idxmin1               <- idx
      }
  }
}

if(plotPoint == T){
    #Create a matrix with all the points-only VMS datapoints
  cc            <- na.omit(as.matrix(cbind(VMS.$declon[pPointsOnly],VMS.$declat[pPointsOnly])))

    #Create smaller matrix to work from for a single interpolation
  timeDiff      <- median(timeDiff,na.rm=T)
  for(poi in 1:dim(cc)[1]){

    iPoiV         <- pPointsOnly[poi]

    res1          <- maxRangeCI(rep(cc[poi,1],2),rep(cc[poi,2],2),timeDiff,rep(VMS.$speed[iPoiV],2))
    
      #First find the boundaries of the mills ellipse, thereafter, add a 10% extra margin, based on the minimum or
      # maximum value. In the longitude direction, take the minimum value, and find the according latitude to go from km to degrees
    res2          <- range(res1[[1]][,1],na.rm=T); boundx <- c(min(res2) - res1[[2]]*0.2*km2Degree(min(res2),
                                                             res1[[1]][which(min(res2)==res1[[1]][,1]),2],1),
                                                             max(res2) + res1[[2]]*0.2*km2Degree(max(res2),
                                                             res1[[1]][which(max(res2)==res1[[1]][,1]),2],1))
    res3          <- range(res1[[1]][,2],na.rm=T); boundy <- c(min(res3) - res1[[2]]*0.2*1/111.2,
                                                             max(res3) + res1[[2]]*0.2*1/111.2)
    if(res1[[3]] == 1){
      boundx <- c(cc[poi,1] - res1[[2]]*0.2*km2Degree(cc[poi,1],cc[poi,2],1),
                  cc[poi,1] + res1[[2]]*0.2*km2Degree(cc[poi,1],cc[poi,2],1))
      boundy <- c(cc[poi,2] - res1[[2]]*0.2*1/111.2,
                  cc[poi,2] + res1[[2]]*0.2*1/111.2)
    }

      #Take the upper right and lower left values as setting the boundings of the smaller matrix
    cc2           <- cbind(boundx,boundy)
    idx           <- na.omit(getGridIndex(cc2,grid,all.inside=T))

      #Work out the other elements of the matrix
    row1          <- min(idx)%/%grid@cells.dim[1]+1;            col1          <- min(idx) - (grid@cells.dim[1]*(row1-1))
    row2          <- max(idx)%/%grid@cells.dim[1]+1;            col2          <- max(idx) - (grid@cells.dim[1]*(row2-1))
    pxheigth      <- abs(row2-row1);                            pxwidth       <- abs(col2-col1)
    bbox          <- matrix(NA,nrow=pxheigth+1,ncol=pxwidth+1); bbox[1,]      <- sort(min(idx) - seq(0,pxwidth,1))
    if(pxheigth > 0) for(i in 1:(pxheigth)) bbox[i+1,] <- bbox[1,] + grid@cells.dim[1]*i
    idx           <- c(bbox)

      #Calculate the distan matrix based on the idx
    distan <- matrix(distance(lon=sPDF@coords[idx,1],lat=sPDF@coords[idx,2],lonRef=cc[poi,1],latRef=cc[poi,2]),
                     nrow=dim(bbox)[1],ncol=dim(bbox)[2])
      #Calculate the distance from begin or endpoint
    linepistan <- distan

      #Reset very small numbers to 0 to get highest values at begin and end point
    linepistan[which(linepistan < 1e-6)]  <- 0
    distan[which(distan < 1e-6)]          <- 0; zeroDistan <- which(distan==0)

    CI                              <- c(matrix(N1p0(distan*params$distscale,0,linepistan^params$sigline,0),ncol=dim(distan)[2],nrow=dim(distan)[1]))
    if(max(CI,na.rm=T) < 0.1) warning("Prediction max(tmpnew) is very small")
    if(length(zeroDistan)>0)  CI[zeroDistan]  <- pmax(CI[zeroDistan],1,na.rm=T)

      #Add point data to final data
    sPDF@data[idx,1]  <- CI + sPDF@data[idx,1]
  }
}
  #Turn point data frame in format that is able to be plotted
sPDF2          <- as(sP,"SpatialPixelsDataFrame")
sPDF2@data     <- data.frame(sPDF@data[,1])
sGDF           <- as(sPDF2,"SpatialGridDataFrame")

if(plot==T){
  library(fields)
  color <- rgb(rep(255,(255*2)),c(rep(255,255),seq(255,0,length.out=255)),c(seq(255,0,length.out=255),rep(0,255)),max=255)
  layout(matrix(c(1,1,1,1,2),nrow=5))
  image(sGDF,col=color)
  mtext("Longitude",side=1,line=3,font=2)
  mtext("Latitude",side=2,line=3,font=2)
  box(); axis(1); axis(2)
  plot(1,1,col="white",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
  legend("center",legend=pretty(seq(0,max(sGDF@data[,1]),length.out=4),n=3),pt.bg=
          color[seq(1,length(color),length.out=length(pretty(seq(0,max(sGDF@data[,1]),length.out=4),n=3)))],
          pch=c(22),col="black",pt.cex=c(2),box.lty=0,ncol=length(pretty(seq(0,max(sGDF@data[,1]),length.out=4),n=3)))

} else { return(sGDF)}
}

