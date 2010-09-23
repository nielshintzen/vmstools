checkSpeedVMS <- function(tacsat=tacsat){

tacsat$SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S")
                          

for(i in 1:(dim(tacsat)[1]-1)){
  if(i==1) VE_REF <- tacsat$VE_REF[1]
  if(VE_REF == tacsat$VE_REF[i+1]){
    dstn  <- degree2Km(tacsat$SI_LONG[i],tacsat$SI_LATI[i],distance(tacsat$SI_LONG[i+1],tacsat$SI_LATI[i+1],tacsat$SI_LONG[i],tacsat$SI_LATI[i]))
    time  <- an(difftime(tacsat$SI_DATIM[i+1],tacsat$SI_DATIM[i],units=c("mins")))
    spd   <- mean(tacsat$SI_SP[i+1],tacsat$SI_SP[i],na.rm=T)
    
    maxdist <- spd*1.852*time/120
    if(dstn > maxdist){
      warning("Speed to low to get at other point")
      tacsat$SI_SP[i] <- NA
    }
  } else {
      VE_REF <- tacsat$VE_REF[i+1]
    }
}
return(tacsat)}
    