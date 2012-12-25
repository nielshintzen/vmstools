intervalTacsat <- function(tacsat,level="trip",weight=c(1,0),fill.na=FALSE){

                    if(length(weight) != 2) stop("weight must be specified as a length 2 numeric vector")
                    weight <- weight / sum(weight,na.rm=TRUE)

                    #- First sort the tacsat dataset
                    tacsat <- sortTacsat(tacsat)

                    #- Add date-time stamp
                    if(!"SI_DATIM" %in% colnames(tacsat)) tacsat$SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")

                    #- If a trip level is specified, the interval rate can be calculated by trip (to make sure that no long
                    #   interval rates between trips occur in comparison to by level="vessel"
                    if(level=="trip"){
                      if(is.null(tacsat$FT_REF)==TRUE) stop("no tripnumber available to merge on trip level")
                      sptacsat      <- split(tacsat,tacsat$VE_REF)
                      tacsat$INTV   <- unlist(lapply(sptacsat,function(x){

                        FT_REF  <- as.factor(x$FT_REF);
                        res     <- by(x,FT_REF,
                                    function(y){
                                      if(nrow(y)>1){
                                        difftime_xmin1  <- c(NA,difftime(y$SI_DATIM[2:nrow(y)],y$SI_DATIM[1:(nrow(y)-1)],units="mins"))
                                        difftime_xplus1 <- c(difftime_xmin1[-1],NA)
                                        if(any(weight == 0)){
                                          if(weight[1] == 0) INTV       <- difftime_xplus1
                                          if(weight[2] == 0) INTV       <- difftime_xmin1
                                        } else {             INTV       <- rowSums(cbind(difftime_xmin1*weight[1],difftime_xplus1*weight[2]))
                                          }
                                        #- If INTV equals NA, then check if there are other possibilities to calculate the interval rate based on a different
                                        #   weight setting.
                                        if(fill.na==TRUE){
                                          idx                           <- which(is.na(INTV)==TRUE)
                                          INTV[idx]                     <- rowSums(cbind(difftime_xmin1[idx],difftime_xplus1[idx]),na.rm=TRUE)
                                          INTV[idx][which(INTV[idx]==0)]<- NA
                                        }

                                        return(INTV)
                                      } else {
                                          return(NA)
                                        }
                                    })
                                                            return(unsplit(res,FT_REF))}))
                      tacsat$INTV[which(tacsat$FT_REF == "0")] <- NA
                    }
                    #- If no trip level is specified, the other option is to calculate interval rate by vessel
                    if(level=="vessel"){
                      difftime_xmin1                  <- c(NA,difftime(tacsat$SI_DATIM[2:nrow(tacsat)],tacsat$SI_DATIM[1:(nrow(tacsat)-1)],units="mins"))
                      difftime_xplus1                 <- c(difftime_xmin1[-1],NA)
                      if(any(weight == 0)){
                        if(weight[1] == 0) INTV       <- difftime_xplus1
                        if(weight[2] == 0) INTV       <- difftime_xmin1
                      } else {             INTV       <- rowSums(cbind(difftime_xmin1*weight[1],difftime_xplus1*weight[2]))
                        }
                      #- If INTV equals NA, then check if there are other possibilities to calculate the interval rate based on a different
                      #   weight setting.
                      if(fill.na==TRUE){
                        idx                           <- which(is.na(INTV)==TRUE)
                        INTV[idx]                     <- rowSums(cbind(difftime_xmin1[idx],difftime_xplus1[idx]),na.rm=TRUE)
                        INTV[idx][which(INTV[idx]==0)]<- NA
                      }
                      tacsat$INTV                     <- INTV

                      vessels                         <- unique(tacsat$VE_REF)
                      first.vessels                   <- unlist(lapply(as.list(vessels),function(x){which(tacsat$VE_REF==x)[1]}))
                      last.vessels                    <- unlist(lapply(as.list(vessels),function(x){rev(which(tacsat$VE_REF==x))[1]}))
                      if(weight[1] != 0) tacsat$INTV[first.vessels]      <- NA
                      if(weight[2] != 0) tacsat$INTV[last.vessels]       <- NA
                      if(fill.na==TRUE) tacsat$INTV[first.vessels]  <- difftime_xplus1[first.vessels]
                      if(fill.na==TRUE) tacsat$INTV[last.vessels]   <- difftime_xmin1[last.vessels]
                    }
                  return(tacsat)}