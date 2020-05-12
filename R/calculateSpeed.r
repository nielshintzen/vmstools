calculateSpeed <- function(tacsat,level="vessel",weight=c(1,1),fill.na=FALSE){

  if(length(weight) != 2) stop("weight must be specified as a length 2 numeric vector")
  weight <- weight / sum(weight,na.rm=TRUE)

  #- Add interval between succeeding points
  tacsatp <- intervalTacsat(tacsat,level=level,weight=weight,fill.na=fill.na)

  if(!"SI_DATIM" %in% colnames(tacsatp)) tacsatp$SI_DATIM  <- as.POSIXct(paste(tacsatp$SI_DATE,  tacsatp$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
  if(level=="trip"){
    if(is.null(tacsatp$FT_REF)==TRUE) stop("no tripnumber available to merge on trip level")
    sptacsat          <- split(tacsatp,tacsatp$VE_REF)
    tacsatp$SI_SPCA   <- unlist(lapply(sptacsat,function(x){

      FT_REF  <- as.factor(x$FT_REF);
      res     <- by(x,FT_REF,
                  function(y){
                    if(nrow(y)>1){
                      interval        <- y$INTV
                      distance_xmin1  <- c(NA,distance(y$SI_LONG[2:nrow(y)],    y$SI_LATI[2:nrow(y)],
                                                       y$SI_LONG[1:(nrow(y)-1)],y$SI_LATI[1:(nrow(y)-1)]))
                      distance_xplus1 <- c(distance_xmin1[-1],NA)
                      if(any(weight == 0)){
                        if(weight[1] == 0) SI_SPCA       <- distance_xplus1 / (interval/60)
                        if(weight[2] == 0) SI_SPCA       <- distance_xmin1 / (interval/60)
                      } else {
                          difftime_xmin1                 <- intervalTacsat(y,level="trip",weight=c(1,0),fill.na=fill.na)$INTV
                          difftime_xplus1                <- intervalTacsat(y,level="trip",weight=c(0,1),fill.na=fill.na)$INTV
                          SI_SPCA                        <- 2*(2* ((distance_xmin1 * weight[1]) / (difftime_xmin1/60)) * (distance_xplus1 * weight[2] / (difftime_xplus1/60)) /
                                                                 (((distance_xmin1 * weight[1]) / (difftime_xmin1/60)) + (distance_xplus1 * weight[2] / (difftime_xplus1/60))))
                        }
                      #- If INTV equals NA, then check if there are other possibilities to calculate the interval rate based on a different
                      #   weight setting.
                      if(fill.na==TRUE){
                        idx                               <- which(is.na(SI_SPCA)==TRUE)
                        if(weight[1] == 0 | weight[2] == 0){
                          speeds                          <- cbind((distance_xmin1[idx]   * weight[1]) / (interval[idx]/60),
                                                                   (distance_xplus1[idx]  * weight[2]) / (interval[idx]/60))
                        } else {
                            speeds                        <- cbind((distance_xmin1[idx]   * weight[1]) / (difftime_xmin1[idx] /60),
                                                                   (distance_xplus1[idx]  * weight[2]) / (difftime_xplus1[idx]/60))
                          }
                        #- Use of 0.5 if ifelse = NA, as counter to 2*...
                        SI_SPCA[idx]                      <- apply(speeds,1,function(x){isfin <- is.finite(x); return(ifelse(all(isfin),mean(x),ifelse(any(isfin),x[isfin],NA)))})
                        #SI_SPCA[idx]                      <- 2*(2* ifelse(is.na(speeds[,1])==TRUE,0.5,speeds[,1]) * ifelse(is.na(speeds[,2])==TRUE,0.5,speeds[,2]) /
                        #                                          (ifelse(is.na(speeds[,1])==TRUE,0,speeds[,1])   + ifelse(is.na(speeds[,2])==TRUE,0,speeds[,2])))
                        if(length(which(SI_SPCA[idx]==0))>0) SI_SPCA[idx][which(SI_SPCA[idx]==0)]<- NA
                      }

                      return(SI_SPCA)
                    } else {
                        return(NA)
                      }
                  })
                                          return(unsplit(res,FT_REF))}))
    tacsatp$SI_SPCA[which(tacsatp$FT_REF == "0")] <- NA

  }
  if(level=="vessel"){
    #- Take the interval for weights c(0,1) and c(1,0) direct from tacsat
    interval                            <- tacsatp$INTV
    distance_xmin1                      <- c(NA,distance(tacsatp$SI_LONG[2:nrow(tacsatp)],     tacsatp$SI_LATI[2:nrow(tacsatp)],
                                                         tacsatp$SI_LONG[1:(nrow(tacsatp)-1)], tacsatp$SI_LATI[1:(nrow(tacsatp)-1)]))
    distance_xplus1                     <- c(distance_xmin1[-1],NA)
    #- Recalculate the interval rate for either outer averages
    difftime_xmin1                      <- intervalTacsat(tacsatp,level="vessel",weight=c(1,0),fill.na=fill.na)$INTV
    difftime_xplus1                     <- intervalTacsat(tacsatp,level="vessel",weight=c(0,1),fill.na=fill.na)$INTV
    if(any(weight == 0)){
      if(weight[1] == 0) SI_SPCA        <- distance_xplus1 / (interval/60)
      if(weight[2] == 0) SI_SPCA        <- distance_xmin1 / (interval/60)
    } else {
        SI_SPCA                         <- 2*(2* ((distance_xmin1 * weight[1]) / (difftime_xmin1/60)) * (distance_xplus1 * weight[2] / (difftime_xplus1/60)) /
                                                (((distance_xmin1 * weight[1]) / (difftime_xmin1/60)) + (distance_xplus1 * weight[2] / (difftime_xplus1/60))))
      }
    #- If INTV equals NA, then check if there are other possibilities to calculate the interval rate based on a different
    #   weight setting.
    if(fill.na==TRUE){
      idx                               <- which(is.na(SI_SPCA)==TRUE)
      if(weight[1] == 0 | weight[2] == 0){
        speeds                          <- cbind((distance_xmin1[idx]   * weight[1]) / (interval[idx]/60),
                                                 (distance_xplus1[idx]  * weight[2]) / (interval[idx]/60))
      } else {
          speeds                        <- cbind((distance_xmin1[idx]   * weight[1]) / (difftime_xmin1[idx] /60),
                                                 (distance_xplus1[idx]  * weight[2]) / (difftime_xplus1[idx]/60))
        }
      #- Use of 0.5 if ifelse = NA, as counter to 2*...
      SI_SPCA[idx]                      <- apply(speeds,1,function(x){isfin <- is.finite(x); return(ifelse(all(isfin),mean(x),ifelse(any(isfin),x[isfin],NA)))})
      #SI_SPCA[idx]                      <- 2*(2* ifelse(is.na(speeds[,1])==TRUE,0.5,speeds[,1]) * ifelse(is.na(speeds[,2])==TRUE,0.5,speeds[,2]) /
      #                                          (ifelse(is.na(speeds[,1])==TRUE,0,speeds[,1])   + ifelse(is.na(speeds[,2])==TRUE,0,speeds[,2]))))
      if(length(which(SI_SPCA[idx]==0))>0) SI_SPCA[idx][which(SI_SPCA[idx]==0)]<- NA
    }
    tacsatp$SI_SPCA                 <- SI_SPCA

    vessels                         <- unique(tacsatp$VE_REF)
    first.vessels                   <- unlist(lapply(as.list(vessels),function(x){which(tacsatp$VE_REF==x)[1]}))
    last.vessels                    <- unlist(lapply(as.list(vessels),function(x){rev(which(tacsatp$VE_REF==x))[1]}))
    if(weight[1] != 0) tacsatp$SI_SPCA[first.vessels]      <- NA
    if(weight[2] != 0) tacsatp$SI_SPCA[last.vessels]       <- NA
    if(fill.na==TRUE) tacsatp$SI_SPCA[first.vessels]          <- distance_xplus1[first.vessels] / (interval[first.vessels]/60)
    if(fill.na==TRUE) tacsatp$SI_SPCA[last.vessels]           <- distance_xmin1[last.vessels] / (interval[last.vessels]/60)

  }
  idx <- which(colnames(tacsatp) %in% c("INTV"))
  #- Convert from km/h to knots / nautical miles an hour
  tacsatp$SI_SPCA     <- tacsatp$SI_SPCA / 1.852
return(if(length(idx)>0){tacsatp[,-idx]}else{tacsatp})}