`findEndVMS` <-
function(VMS
                          ,startVMS #Starting point of VMS
                          ,interval #Specify in minutes, NULL means use all points  
                          ,margin   #Specify the margin in minutes it might deviate from the interval time, in minutes
                       ){

clStartVMS  <- startVMS #Total VMS list starting point instead of subset use
iShip       <- VMS$ship[startVMS]
VMS.        <- subset(VMS,ship==iShip)
startVMS    <- which(VMS$ship[startVMS] == VMS.$ship & VMS$date[startVMS] == VMS.$date)
if(clStartVMS != dim(VMS)[1]){
  if(VMS$ship[clStartVMS] != VMS$ship[clStartVMS+1]){
      #End of dataset reached
    endDataSet <- 1
    endVMS <- NA
  } else {
        #Calculate the difference in time between the starting VMS point and its succeeding points
      diffTime  <- difftime(VMS.$date[(startVMS+1):dim(VMS.)[1]],VMS.$date[startVMS],units=c("mins"))
      if(length(which(diffTime >= (interval-margin) & diffTime <= (interval+margin)))==0){
        warning("No succeeding point found, no interpolation possible")
        endVMS  <- NA
          #Check if end of dataset has been reached
        ifelse(all((diffTime < (interval-margin))==T),endDataSet <- 1,endDataSet <- 0)
      } else {
          res <- which(diffTime >= (interval-margin) & diffTime <= (interval+margin))
          if(length(res)>1){
            res2        <- which.min(abs(interval-an(diffTime[res])))
            endVMS      <- startVMS + res[res2]
            endDataSet  <- 0
          } else {
              endVMS      <- startVMS + res
              endDataSet  <- 0
            }
          }
        #Build-in check
      if(is.na(endVMS)==F){
        if(!an(difftime(VMS.$date[endVMS],VMS.$date[startVMS],units=c("mins"))) %in% seq((interval-margin),(interval+margin),1)) stop("found endVMS point not within interval range")
        endVMS <- clStartVMS + (endVMS - startVMS)
      }

    }
} else { endDataSet <- 1; endVMS <- NA}

return(c(endVMS,endDataSet))}

