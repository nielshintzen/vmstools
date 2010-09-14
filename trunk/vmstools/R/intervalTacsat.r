intervalTacsat <- function(tacsat){
                    tacsat$DATIM  <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M:%S")
                    sptacsat      <- split(tacsat,tacsat$VE_REF)
                    tacsat$INTV   <- unlist(lapply(sptacsat,function(x){FT_REF  <- as.factor(x$FT_REF);
                                                                        res     <- by(x,FT_REF,function(y){ if(nrow(y)>1){ return(c(NA,difftime(y$DATIM[2:nrow(y)],y$DATIM[1:(nrow(y)-1)],units="mins")))}else{return(NA)}})
                                                            return(unlist(res))}))
                  return(tacsat)}