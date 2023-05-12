#' Calculate the distance (in km) from a Tacsat dataset
#' 
#' Calculates the distance of a Tacsat dataset with specification of succeeding
#' datapoints. Distance is only calculated between these specified datapoints
#' 
#' index is designed as a matrix where rows represent succeeding datapoints,
#' column 1 represent start points, column 2 represent end points.
#' 
#' @param tacsat tacsat (normal or high ping rate) dataset
#' @param index Matrix with specification of succeeding datapoints (see details
#' for format)
#' @author Niels T. Hintzen
#' @seealso \code{\link{distance}}, \code{\link{distanceInterpolation}},
#' \code{\link{diffInter}}
#' @references EU lot 2 project
#' @examples
#' 
#' 
#' data(tacsat)
#' #Speed threshold points (two values), NULL means use all points
#' st        <- c(2,6)
#' #Remove duplicate records in VMS dataset
#' remDup    <- TRUE
#' 
#'   #Sort the VMS data
#' tacsat     <- sortTacsat(tacsat)
#' tacsat     <- tacsat[1:1000,]
#' 
#'   #Filter the VMS data
#' tacsat     <- filterTacsat(tacsat,st,NULL,remDup)
#' 
#' distanceTacsat(tacsat,matrix(c(2,3,3,4),nrow=2,ncol=2,
#'                 dimnames=list(1:2,c("startpoint","endpoint"))))
#' #6.335944 14.847291
#' 
#' @export distanceTacsat
distanceTacsat <- function(tacsat,index){

                   res <- unlist(lapply(as.list(1:dim(index)[1]),function(x){
                              iS  <- index[x,1]
                              iE  <- index[x,2]
                              iL  <- iE-iS+1
                              res <- distance(tacsat[iS:iE,]$SI_LONG[2:iL],tacsat[iS:iE,]$SI_LATI[2:iL],
                              tacsat[iS:iE,]$SI_LONG[1:(iL-1)],tacsat[iS:iE,]$SI_LATI[1:(iL-1)])
                           return(sum(res,na.rm=TRUE))}))

               return(res)}
