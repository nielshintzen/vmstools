checkSpeedVMS <- function(VMS){

for(i in 1:(dim(VMS)[1]-1)){
  if(i==1) ship <- VMS$ship[1]
  if(ship == VMS$ship[i+1]){
    dstn  <- degree2Km(VMS$declon[i],VMS$declat[i],distance(VMS$declon[i+1],VMS$declat[i+1],VMS$declon[i],VMS$declat[i]))
    time  <- an(difftime(VMS$date[i+1],VMS$date[i],units=c("mins")))
    spd   <- mean(VMS$speed[i+1],VMS$speed[i],na.rm=T)
    
    maxdist <- spd*1.852*time/120
    if(dstn > maxdist){
      warning("Speed to low to get at other point")
      VMS$speed[i] <- NA
    }
  } else {
      ship <- VMS$ship[i+1]
    }
}
return(VMS)}
    