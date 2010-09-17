intervalTacsat <- function(tacsat,level="trip"){
                    tacsat$DATIM  <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M:%S")
                    if(level=="trip"){
                      if(is.null(tacsat$FT_REF)==T) stop("no tripnumber available to merge on trip level")
                      sptacsat      <- split(tacsat,tacsat$VE_REF)
                      tacsat$INTV   <- unlist(lapply(sptacsat,function(x){FT_REF  <- as.factor(x$FT_REF);
                                                                        res     <- by(x,FT_REF,function(y){ if(nrow(y)>1){ return(c(NA,difftime(y$DATIM[2:nrow(y)],y$DATIM[1:(nrow(y)-1)],units="mins")))}else{return(NA)}})
                                                            return(unlist(res))}))
                    }
                    if(level=="vessel"){
                      tacsat$INTV[2:nrow(tacsat)]  <- difftime(tacsat$DATIM[2:nrow(tacsat)],tacsat$DATIM[1:(nrow(tacsat)-1)],units="mins")
                      vessels                         <- unique(tacsat$VE_REF)
                      first.vessels                   <- unlist(lapply(as.list(vessels),function(x){which(tacsat$VE_REF==x)[1]}))
                      tacsat$INTV[first.vessels-1] <- NA
                    }
                  return(tacsat)}