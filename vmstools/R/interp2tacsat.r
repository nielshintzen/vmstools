
interp2tacsat <- function(interpolation,tacsat,npoints=10){

# This function takes the list of tracks output by interpolateTacsat and converts those back to tacsat format.

interpolationEQ <- equalDistance(interpolation,npoints)  #Divide points equally along interpolated track (default is 10).

res <- lapply(interpolationEQ,function(x){
                                  idx                     <- x[1,]; x <- data.frame(x)
                                  colnames(x)             <- c("SI_LONG","SI_LATI")
                                  x$VE_COU                <- rep(tacsat$VE_COU[idx[1]],nrow(x))
                                  x$VE_REF                <- rep(tacsat$VE_REF[idx[1]],nrow(x))
                                  x$FT_REF                <- rep(tacsat$FT_REF[idx[1]],nrow(x))
                                  x$SI_DATE               <- tacsat$SI_DATE[idx[1]]
                                  x$SI_TIME               <- tacsat$SI_TIME[idx[1]]
                                  x$SI_DATIM              <- tacsat$SI_DATIM[idx[1]]
                                  x$SI_DATIM[-c(1:2)]     <- as.POSIXct(cumsum(rep(difftime(tacsat$SI_DATIM[idx[2]],tacsat$SI_DATIM[idx[1]],unit="secs")/(nrow(x)-2),nrow(x)-2))+tacsat$SI_DATIM[idx[1]],tz="GMT",format = "%d/%m/%Y  %H:%M")
                                  x$SI_SP                 <- mean(tacsat$SI_SP[idx[1]],tacsat$SI_SP[idx[2]],na.rm=T)
                                  x$SI_HE                 <- NA;
                                  x$SI_HE[-c(1,nrow(x))]  <- bearing(x$SI_LONG[3:nrow(x)],x$SI_LATI[3:nrow(x)],x$SI_LONG[2:(nrow(x)-1)],x$SI_LATI[2:(nrow(x)-1)])
                                return(x[-c(1,2,nrow(x)),])})

interpolationTot  <- do.call(rbind,res)
tacsatInt         <- rbind(interpolationTot,tacsat[,colnames(interpolationTot)])
tacsatInt         <- sortTacsat(tacsatInt)

return(tacsatInt)

}
