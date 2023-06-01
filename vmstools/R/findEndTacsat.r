#' Finding the succeeding Tacsat datapoint based on an interval with a
#' specified margin
#' 
#' To create an interpolation, two succeeding Tacsat datapoints are needed.
#' This function finds the succeeding Tacsat point and tests if the point is
#' within the specified time interval and margins. As well, if no succeeding
#' datapoint can be found, this information is returned
#' 
#' Interval: In most Tacsat datasets the succeeding datapoint can be found 1 or
#' 2 hours appart. This interval time should be specified here. Interval can
#' also be specified as e.g. 15 minutes if the Tacsat / GPS dataset allows
#' this. Margin: Hardly ever, the interval time is precise. To allow some
#' deviation from the 1 or 2 hour interval the margin can be adjusted.
#' 
#' The result returned consists of 2 values. The first value is the index of
#' the Tacsat set specified of the succeeding datapoint. The second value
#' indicates if the dataset has ended. If 1st: NA and 2nd 0 then no succeeding
#' Tacsat point could be found in the specified interval. If 1st: NA and 2nd 1
#' then no succeeding Tacsat point could be found and end of dataset for a
#' specific vessel has been reached. If 1st: NA and 2nd 2 then no succeeding
#' Tacsat point could be found and end of complete dataset has been reached. If
#' 1st: value then 2nd will be 0, succeeding Tacsat point is found and is
#' specified in 1st value.
#' 
#' @param tacsat The Tacsat dataset
#' @param startTacsat Index of Tacsat dataset of startpoint of interpolation
#' @param interval Time in minutes between the succeeding datapoints
#' @param margin Deviation from specified interval to find succeeding
#' datapoints
#' @note This function is called inside interpolateTacsat()
#' @author Niels T. Hintzen
#' @seealso \code{\link{filterTacsat}}, \code{\link{interpolateTacsat}}
#' @references EU lot 2 project
#' @examples
#' 
#' data(tacsat)
#' startTacsat <- 2
#' #result: 3 0 Succeeding point = tacsat[3,]
#' # and end dataset has not been reached yet.
#' findEndTacsat(tacsat,startTacsat,interval=120,margin=10)
#' 
#' @export findEndTacsat
`findEndTacsat` <-
function(SI_DATIM
                          ,startVMS #Starting point of VMS
                          ,interval #Specify in minutes, NULL means use all points  
                          ,margin   #Specify the margin in minutes it might deviate from the interval time, in minutes
                       ){

  #Calculate the difference in time between the starting VMS point and its succeeding points
  diffTime          <- difftime(SI_DATIM[(startVMS+1):length(SI_DATIM)],SI_DATIM[startVMS],units=c("mins"))
  if(length(which(diffTime >= (interval-margin) & diffTime <= (interval+margin)))==0){
    warning("No succeeding point found, no interpolation possible")
    endVMS          <- NA
    endDataSet      <- 3
      #Check if end of dataset has been reached
    ifelse(all((diffTime < (interval-margin))==TRUE),endDataSet <- 1,endDataSet <- 0)
  } else {
      res           <- which(diffTime >= (interval-margin) & diffTime <= (interval+margin))
      if(length(res)>1){
        res2        <- which.min(abs(interval-an(diffTime[res])))
        endVMS      <- startVMS + res[res2]
        endDataSet  <- 0
      } else {
          endVMS    <- startVMS + res
          endDataSet<- 0
        }
      }
    #Build-in check
  if(is.na(endVMS)==FALSE){
    if(!an(difftime(SI_DATIM[endVMS],SI_DATIM[startVMS],units=c("mins"))) %in% seq((interval-margin),(interval+margin),1)) stop("found endVMS point not within interval range")
    endVMS <- (endVMS - startVMS)
  }
return(c(endVMS,endDataSet))}

