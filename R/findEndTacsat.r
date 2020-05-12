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

