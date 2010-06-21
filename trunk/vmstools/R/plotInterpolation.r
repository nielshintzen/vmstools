`plotInterpolation` <-
function(interpolation
                                           ,VMS
                                           ,grid
                                           ,plot=F
                                           ,singlePoints=F
                                           ){

VMS. <- VMS
                                           
  #Define the single points, and the interpolated points
pPoints     <- seq(1,dim(VMS.)[1],1)
pInterOnly  <- unique(unlist(lapply(interpolation,function(x){return(x[1,])})))
pPointsOnly <- which(is.na(pmatch(pPoints,pInterOnly))==T)

  #Create SpatialGridDataFrame, the object where we can store the impact data
spatialGrid <- SpatialGrid(grid=grid)
gridded(spatialGrid) = TRUE
sP          <- as(spatialGrid,"SpatialPixels")
sPDF        <- as(sP,"SpatialPixelsDataFrame")
sPDF@data   <- data.frame(rep(0,length(sPDF@grid.index)))

if(singlePoints == T){
    #Create a matrix with all the points-only VMS datapoints
  cc            <- na.omit(as.matrix(cbind(VMS.$declon[pPointsOnly],VMS.$declat[pPointsOnly])))
  
    #Overlay the points to the grid
  idx               <- na.omit(getGridIndex(cc,grid,all.inside=T))
  sPDF@data[an(dimnames(table(idx))$idx),1]  <- an(table(idx))
}
  #Overlay the interpolations to the grid
res1          <- matrix(unlist(interpolation),nrow=dim(interpolation[[1]])[1],ncol=2*length(interpolation),
                        dimnames=list(c("Point",seq(1,dim(interpolation[[1]])[1]-1,1)),rep(c("x","y"),length(interpolation))))
uniquePoints  <- which(duplicated(unlist(lapply(interpolation,function(x){return(x[1,])})))==F)
res5          <- sort(c(cbind(uniquePoints[which(uniquePoints%%2==1)],uniquePoints[which(uniquePoints%%2==1)]+1)))
res6          <- which(!seq(1,dim(res1)[2],1) %in% res5)
res7          <- cbind( c(c(res1[-1,res5][,which(dimnames(res1[-1,res5])[[2]]=="x")]),c(res1[-c(1,2),res6][,which(dimnames(res1[-c(1,2),res6])[[2]]=="x")])),
                        c(c(res1[-1,res5][,which(dimnames(res1[-1,res5])[[2]]=="y")]),c(res1[-c(1,2),res6][,which(dimnames(res1[-c(1,2),res6])[[2]]=="y")])))
  #Create a matrix with all the interpolation-only VMS datapoints
cc            <- na.omit(res7)
idx           <- na.omit(getGridIndex(cc,grid,all.inside=T))

idxagg        <- table(idx)
idxagg2       <- cbind(an(dimnames(idxagg)[[1]]),idxagg)
sPDF@data[idxagg2[,1],1] <- sPDF@data[idxagg2[,1],1] + idxagg2[,2]

  #Turn point data frame in format that is able to be plotted
sGDF              <- as(sPDF,"SpatialGridDataFrame")
  #Plot the results
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

