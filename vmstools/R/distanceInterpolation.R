#' Calculate the distance (in km) from an interpolated dataset
#' 
#' For each interpolation between two succeeding datapoints, the distance
#' travelled is computed based on the function 'distance()'. This function is
#' an easy to use wrapper for a whole interpolated dataset
#' 
#' On default, each interpolation consists of 100 points. Hence, 99 distances
#' are computed for each interpolation. Returned are 5
#' 
#' @param interpolation Interpolated dataset as output from the function
#' 'interpolateVMS'
#' @author Niels T. Hintzen
#' @seealso \code{\link{distance}}, \code{\link{distanceTacsat}}
#' @references EU lot 2 project
#' @examples
#' 
#' data(tacsat)
#' #Speed threshold points (two values), NULL means use all points
#' st        <- c(2,6)
#' #Remove duplicate records in VMS dataset
#' remDup    <- TRUE
#' 
#'   #Sort the tacsat data
#' tacsat     <- sortTacsat(tacsat)
#' tacsat     <- tacsat[1:1000,]
#' 
#'   #Filter the tacsat data
#' tacsat     <- filterTacsat(tacsat,st,NULL,remDup)
#' 
#' interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
#'                   res=100,method="cHs",params=list(fm=0.3,distscale=20,
#'                   sigline=0.2,st=c(2,6)),headingAdjustment=0)
#'   #Number of interpolations:
#' length(interpolation)
#'   #Calculate distance:
#'   #Returns 159 values, which represents the distance travelled in
#'   #each of the 159 interpolations
#' distanceInterpolation(interpolation)
#' 
#' @export distanceInterpolation
distanceInterpolation <- function(interpolation){
                                               
                            res <- unlist(lapply(interpolation,function(x){
                                                 dims        <- dim(x)
                                                 res         <- distance(x[3:dims[1],1],x[3:dims[1],2],x[2:(dims[1]-1),1],x[2:(dims[1]-1),2])
                                                 return(sum(res,na.rm=TRUE))}))
                     
                            return(res)}


  


  
