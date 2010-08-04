outsideMaxRange <- function(intLon
                                  ,intLat
                                  ,vmsIdx1
                                  ,vmsIdx2
                                  ,VMS.
                                  ,grid
                                  ,sPDF
                                  ,interpolation
                                  ,int
                                  ,params){
                                  
  #Calculate the confidence interval                               
res <- calculateCI(intLon,intLat,vmsIdx1,vmsIdx2,VMS.,grid,sPDF,interpolation,int,params)
CI  <- res[[1]] #computed CI
idx <- res[[2]] #identifier of spatial data frame index
mxR <- res[[3]] #max range
grid<- res[[4]] #new defined grid if expanded
sPDF<- res[[5]] #new defined grid if expanded

  #Write the confidence interval in the spatial data frame structure
sPDF@data[idx,1]  <- CI

#sPDF2          <- as(sP,"SpatialPixelsDataFrame")
#sPDF2@data     <- data.frame(sPDF@data[,1])
#sGDF           <- as(sPDF2,"SpatialGridDataFrame")
#
#image(sGDF)
#lines(mxR[[1]][,1],mxR[[1]][,2])
#points(sPDF@coords[idx,1],sPDF@coords[idx,2],pch=19,cex=0.3,col="red")
#points(sPDF@coords[idx,1][which(res2 == 0)],sPDF@coords[idx,2][which(res2 == 0)],pch=19,cex=0.3,col="blue")
#

  #Calculate which proportions of the CI are located outside and inside the maximum range ellipse
  #As well, calculate the maximum value of the CI as these should represent 1 at both the begin and end VMS data points
res2      <- point.in.polygon(sPDF@coords[idx,1],sPDF@coords[idx,2],mxR[[1]][,1],mxR[[1]][,2])
insideR   <- sum(sPDF@data[idx,1][which(res2 == 1)],na.rm=T) #Sum of total CI values inside the maximum range, should ideally be all the grid cells with values
outsideR  <- sum(sPDF@data[idx,1][which(res2 == 0)],na.rm=T) #Sum of total CI values outside the maximum range, should ideally be 0
maxR      <- max(sPDF@data[idx,1][which(res2 == 1)],na.rm=T) #Top of the CI, should ideally equal to 1

return(list(insideR,outsideR,maxR))}
