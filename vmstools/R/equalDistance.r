#' Interpolated points at equal distance
#' 
#' Returns the interpolated dataset with only those points remaining that are
#' at equal eucledian distance from each other, with the number of points to
#' retreive remaining.
#' 
#' 
#' @param interpolation interpolated dataset obtained from the interpolation()
#' function
#' @param res number of points to retreive from function
#' @author Niels T. Hintzen
#' @seealso \code{\link{filterTacsat}}, \code{\link{tacsat}},
#' \code{\link{interpolateTacsat}}
#' @references EU lot 2 project
#' @examples
#' 
#' data(tacsat)
#' 
#' #Sort the VMS data
#' tacsat     <- sortTacsat(tacsat)
#' tacsat     <- tacsat[1:1000,]
#' 
#' #Filter the Tacsat data
#' tacsat     <- filterTacsat(tacsat,st=c(2,6),hd=NULL)
#' 
#' #Interpolate the VMS data
#' interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
#'                     res=100,method="cHs",params=list(fm=0.5,distscale=20,
#'                     sigline=0.2,st=c(2,6)),headingAdjustment=0)
#' 
#' #Get a set back with only 10 points per interpolation at equal distance
#' ed_interpolation <- equalDistance(interpolation,10)
#' 
#' @export equalDistance
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
