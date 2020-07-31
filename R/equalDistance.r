equalDistance <- function(interpolation,res=10){

                    #Calculate ditance of all interpolations at the same time
                    totDist <- distanceInterpolation(interpolation)
                    #Get dimensions of interpolations
                    lngInt  <- lapply(interpolation,dim)

                    #Warn if resolution of equal distance is too high compared to original resolution of interpolation
                    if(min(unlist(lngInt)[seq(1,length(totDist),2)],na.rm=TRUE) < 9*res) warnings("Number of intermediate points in the interpolation might be too small for the equal distance pionts chosen")

                    #Get distance steps to get equal distance
                    eqStep  <- totDist/(res-1)
                    
                    #Get x-y values of all interpolations
                    intidx  <- matrix(unlist(lapply(interpolation,function(x){return(x[1,])})),ncol=2,byrow=TRUE)

                    #Do the calculation
                    result  <- lapply(interpolation,function(ind){
                                                      i       <- which(intidx[,1] == ind[1,1] & intidx[,2] == ind[1,2])
                                                      idx     <- apply(abs(outer(
                                                                                 cumsum(distance(ind[3:lngInt[[i]][1],1],ind[3:lngInt[[i]][1],2],ind[2:(lngInt[[i]][1]-1),1],ind[2:(lngInt[[i]][1]-1),2])),
                                                                                 seq(eqStep[i],totDist[i],eqStep[i]),
                                                                                 "-")),
                                                                       2,which.min)+1
                                                      idx     <- c(1,idx)
                                                    return(ind[c(1,idx+1),])})
                  #Return the equal distance interpolated set in the same format as the interpolated dataset (as a list)
                  return(result)}
