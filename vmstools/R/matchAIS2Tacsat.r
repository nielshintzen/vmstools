matchAIS2Tacsat <- function(ais,tacsat,interval=120,margin=12,unitAIS="mins"){

  #- Match vessels
  selectVessels <- unique(tacsat$VE_REF)[which(unique(tacsat$VE_REF) %in% unique(ais$VE_REF))]


  if(class(tacsat$SI_DATIM)[1] != "POSIXct")
    tacsat$SI_DATIM  <- as.POSIXct(tacsat$SI_DATIM,format="%Y-%m-%d %H:%M:%S")
  if(class(ais$SI_DATIM)[1] != "POSIXct")
    ais$SI_DATIM     <- as.POSIXct(ais$SI_DATIM,   format="%Y-%m-%d %H:%M:%S")

  #- Generate ID
  ais$ID        <- 1:nrow(ais)
  tacsat$ID     <- 1:nrow(tacsat)
  
  #- Make subset to match vessels
  subais        <- subset(ais,VE_REF %in% selectVessels)
  subtacsat     <- subset(tacsat,VE_REF %in% selectVessels)

  #- Match days
  subtacsat     <- sortTacsat(subtacsat)
  subais        <- sortTacsat(subais)

  if(class(subtacsat$SI_DATIM)[1] != "POSIXct")
    subtacsat$SI_DATIM  <- as.POSIXct(subtacsat$SI_DATIM,format="%Y-%m-%d %H:%M:%S")
  if(class(subais$SI_DATIM)[1] != "POSIXct")
    subais$SI_DATIM     <- as.POSIXct(subais$SI_DATIM,   format="%Y-%m-%d %H:%M:%S")
  #- Round timestamp
  subtacsat$SI_DATIM <- round(subtacsat$SI_DATIM,unitAIS)
  subais$SI_DATIM    <- round(subais$SI_DATIM,unitAIS)

  subtacsat$jul <- julian(as.Date(subtacsat$SI_DATIM))
  subais$jul    <- julian(as.Date(subais$SI_DATIM))
  matchJul      <- unique(subtacsat$jul)[which(unique(subtacsat$jul) %in% unique(subais$jul))]

  subtacsat     <- subset(subtacsat,jul %in% matchJul)
  subais        <- subset(subais,jul %in% matchJul)

  #- Split by vessel
  spltTacsat    <- split(subtacsat,subtacsat$VE_REF)
  spltAis       <- split(subais,subais$VE_REF)
  
  #- Link AIS to VMS by vessel
  matchAIS      <- list()
  matchTAC      <- list()
  for(iRef in names(spltTacsat)){
    vesAIS      <- spltAis[[iRef]]
    vesTAC      <- spltTacsat[[iRef]]
    vesTAC      <- intervalTacsat(vesTAC,level="vessel",fill.na=T)
    datimAIS    <- vesAIS$SI_DATIM

    #- Get time-frames in between VMS is valid
    idx         <- which(vesTAC$INTV >= (interval - margin) & vesTAC$INTV <= (interval + margin))
    if(length(idx)>0){
      if(nrow(vesTAC) %in% idx)
        idx       <- idx[-length(idx)]
      pairsTAC    <- as.data.frame(cbind(an(vesTAC[idx,"SI_DATIM"]),an(vesTAC[idx+1,"SI_DATIM"])))
      pairsTAC$FT_INTV <- 1:nrow(pairsTAC)
      #pairsTAC1    <- cbind.data.frame(vesTAC[idx,"SI_DATIM"],vesTAC[idx+1,"SI_DATIM"])
      if(unitAIS=="mins")
        res         <- do.call(rbind,apply(pairsTAC,1,FUN=function(x){return(cbind(seq(x[1],x[2],by=60),x[3]))}))
      if(unitAIS=="hours")
        res         <- do.call(rbind,apply(pairsTAC,1,FUN=function(x){return(cbind(seq(x[1],x[2],by=60*60),x[3]))}))
      if(unitAIS=="days")
        res         <- do.call(rbind,apply(pairsTAC,1,FUN=function(x){return(cbind(seq(x[1],x[2],by=60*60*24),x[3]))}))

      #- Find AIS pings in between time-frames
      submatch         <- which(an(vesAIS$SI_DATIM) %in% res[,1])
      matchAIS[[iRef]] <- cbind(vesAIS[submatch,],res[match(an(vesAIS[submatch,"SI_DATIM"]),res[,1]),2])
      colnames(matchAIS[[iRef]])[ncol(matchAIS[[iRef]])] <- "FT_INTV"
      matchTAC[[iRef]] <- cbind(vesTAC[idx,],pairsTAC$FT_INTV)
    }
  }
  matchAIS      <- do.call(rbind,matchAIS)
  matchTAC      <- do.call(rbind,matchTAC)
  colnames(matchTAC)[ncol(matchTAC)] <- "FT_INTV"
return(list(ais   =matchAIS,
            tacsat=matchTAC))}


#- Example
#load("./ais.RData")
#load("./tacsat.RData")
##- Make sure the columns INTV and SI_DATIM are in both datasets
#
#output <- matchAIS2Tacsat(ais,tacsat,interval=120,margin=12,unitAIS="mins")
#
##- Write file
#filteredAIS <- output[["ais"]]
#filteredTAC <- output[["tacsat"]]
#effBySeqAIS <- aggregate(filteredAIS$INTV,by=list(filteredAIS$VE_REF,filteredAIS$FT_INTV),FUN=sum,na.rm=T)
#effBySeqTAC <- aggregate(filteredTAC$INTV,by=list(filteredTAC$VE_REF,filteredTAC$FT_INTV),FUN=sum,na.rm=T)
#colnames(effBySeqAIS) <- c("VE_REF","FT_INTV","effAIS")
#colnames(effBySeqTAC) <- c("VE_REF","FT_INTV","effTAC")
#effBySeqComb<- merge(effBySeqAIS,effBySeqTAC)
#effBySeqComb$ratio    <- effBySeqComb$effAIS / effBySeqComb$effTAC
#
##- Take a particular threshold, I've set it very small here but I suggest to
##  find an optimum in which you don't throw out too much data
#idx         <- which(effBySeqComb$ratio > 0.90 & effBySeqComb$ratio < 1.10)
#seq2select  <- apply(effBySeqComb[idx,c("VE_REF","FT_INTV")],1,paste,collapse="_")
#filteredAIS$VE_REF_INTV <- apply(filteredAIS[,c("VE_REF","FT_INTV")],1,paste,collapse="_")
#filteredTAC$VE_REF_INTV <- apply(filteredTAC[,c("VE_REF","FT_INTV")],1,paste,collapse="_")
#aisnew      <- subset(filteredAIS,VE_REF_INTV %in% seq2select)
#tacsatnew   <- subset(filteredTAC,VE_REF_INTV %in% seq2select)
#
##- These two values should be very close together
#sum(aisnew$INTV)
#sum(tacsatnew$INTV)
    
    
    
    
    
    
    
    
    