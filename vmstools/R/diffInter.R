#' Calculate distance between interpolation and reference set
#' 
#' Calculate for each interpolation the distance at fixed points, depending on
#' the number of points present in the reference set, between the interpolated
#' trajectory and the reference trajectory. This indicates the deviation of the
#' interpolated set to the reference set.
#' 
#' Each interpolation has a start and end point which are similar to the start
#' and end point of the reference set. In between these two points, the
#' reference set can have more in-between points which are not used for
#' interpolation, depending on the interval chosen in the 'interpolateVMS'
#' function. These points are matched up with the points on the interpolated
#' track. For each of these reference in-between points, it is calculated which
#' fraction of the total distance was travelled. The same fraction is applied
#' to the interpolated track, and the point that matches with this distance is
#' linked to the reference set point. The distance function is used to
#' calculate the distance between the reference and interpolated dataset.
#' 
#' The mean, logmean, standard deviation, sdlog and total sum of these
#' distances is returned as a matrix.
#' 
#' @param interpolation Interpolated dataset as output from the function
#' 'interpolateVMS'
#' @param reference Original high higher resolution VMS dataset
#' @author Niels T. Hintzen
#' @seealso \code{\link{distance}}, \code{\link{distanceInterpolation}},
#' \code{\link{distanceVMS}}
#' @references Hintzen et al. 2010 Improved estimation of trawling tracks using
#' cubic Hermite spline interpolation of position registration data, EU lot 2
#' project
#' @examples
#' 
#' \dontrun{
#' data(VMShf)
#' 
#' #-Put the data in the right format
#' VMShf$SI_DATE <- format(as.Date(VMShf$date),"%d/%m/%Y")
#' VMShf$SI_TIME <- format(VMShf$date, "%H:%M")
#' 
#' colnames(VMShf) <- c("VE_REF","SI_LATI","SI_LONG","SI_SP",
#'                      "SI_HE","SI_DATIM","SI_DATE","SI_TIME")
#' 
#' #-Sort the data and remove non-fishing pings
#' VMShf   <- sortTacsat(VMShf)
#' VMShf   <- filterTacsat(VMShf,c(2,6),NULL,T)
#' 
#' interpolation <- interpolateTacsat(VMShf,interval=120,margin=10,res=100,
#'                     method="cHs",params=list(fm=0.3,distscale=20,
#'                     sigline=0.2,st=c(2,6)),headingAdjustment=0)
#' 
#' #Returns a matrix with 21 rows (each row represents 1 interpolation)
#' # and 5 measures
#' cHs <- diffInter(interpolation,VMShf)
#' }
#' 
#' @export diffInter
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
