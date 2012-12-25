diffInter <- function(interpolation
                                   ,reference){


      #Get the starting and ending positions of the interpolations
    interIdx    <- matrix(unlist(lapply(interpolation,function(x){return(x[1,])})),ncol=2,nrow=length(interpolation),dimnames=list(interpolation=1:length(interpolation),c("x","y")),byrow=TRUE)
      #Store the deviations from the interpolation and reference dataset
    storeDiffs  <- matrix(NA,nrow=length(interpolation),ncol=5,dimnames=list(1:length(interpolation),c("mean","logmean","sd","logsd","sum")))

      #Loop over all the interpolations
    for(i in 1:length(interpolation)){
      int       <- interpolation[[i]]
      ref       <- reference[seq(interIdx[i,1],interIdx[i,2],1),]

        #Calculate the difference between each datapoint
      distInt <- distance(int[3:dim(int)[1],1],  int[3:dim(int)[1],2],  int[2:(dim(int)[1]-1),1],      int[2:(dim(int)[1]-1),2])
      distRef <- distance(ref$SI_LONG[2:dim(ref)[1]], ref$SI_LATI[2:dim(ref)[1]], ref$SI_LONG[1:(dim(ref)[1]-1)], ref$SI_LATI[1:(dim(ref)[1]-1)])
        #To calculate the total distance travelled, sum all individual distances
      cumsumInt     <- cumsum(distInt)
      cumsumRef     <- cumsum(distRef)
        #To select the points on the interpolated track to match up with the reference points we do the following
        #- First rescale the distance travelled within the reference set to equal the total distance in the interpolated set
        #- Than substract the distance travelled in the reference set for each point from the distance travelled in the interpolated set, between all points
        #- Search for the point in the interpolated set that comes closest to that distance travelled
        # By doing this, we assume that an equal distance of the total in the reference set and in the interpolated set is travelled
        # This enables the user too to make use of irregular reference set polling rates, as matching points are found based on distance travelled
        # As well, if vessels speed between two points and slow down between other, this does match up better with the interpolated set
      matchRefDist  <- c(cumsumRef / (rev(cumsumRef)[1]/rev(cumsumInt)[1]))
      matchRef      <- c(1,apply(abs(outer(matchRefDist,cumsumInt,"-")),1,which.min)+1)

      matchPx       <- int[matchRef+1,1]
      matchPy       <- int[matchRef+1,2]

        #Calculate the distance between the reference points and the points on the interpolated track that are matched
      res <- distance(matchPx,matchPy,ref$SI_LONG,ref$SI_LATI)
        #Store the differences for each interpolation
      storeDiffs[i,]<- c( mean(res[-c(1,length(res))],na.rm=TRUE),
                          exp(mean(log(res[-c(1,length(res))])[which(is.finite(log(res[-c(1,length(res))]))==TRUE)],na.rm=TRUE)),
                          sd(res[-c(1,length(res))],na.rm=TRUE),
                          exp(sd(log(res[-c(1,length(res))])[which(is.finite(log(res[-c(1,length(res))]))==TRUE)],na.rm=TRUE)),
                          sum(res[-c(1,length(res))],na.rm=TRUE))
    }
return(storeDiffs)}
