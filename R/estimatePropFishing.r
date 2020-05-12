estimatePropFishing <- function(tacsat,eflalo,by=c("LE_GEAR","VE_REF")){

  if(!"SI_STATE" %in% colnames(tacsat))
    stop("Provide 'SI_STATE' column in tacsat dataset")


  #-----------------------------------------------------------------------------
  #- Add columns on dates
  #-----------------------------------------------------------------------------
  if(!"SI_DATIM" %in% colnames(tacsat))
    tacsat$SI_DATIM <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M")
  if(!"FT_DDATIM" %in% colnames(eflalo))
    eflalo$FT_DDATIM <- as.POSIXct(paste(eflalo$FT_DDAT,eflalo$FT_DTIME),format="%d/%m/%Y %H:%M",tz = "GMT")
  if(!"FT_LDATIM" %in% colnames(eflalo))
    eflalo$FT_LDATIM <- as.POSIXct(paste(eflalo$FT_LDAT,eflalo$FT_LTIME),format="%d/%m/%Y %H:%M",tz = "GMT")
  if(!"INTV" %in% colnames(tacsat))
    stop("Specify time interval column in tacsat (e.g. use intervalTacsat)")

  #- Add date notation
  if("SI_DAY" %in% unique(c(by))){
    eflalo$SI_DAY <- yday(as.POSIXct(eflalo$LE_CDAT,format="%d/%m/%Y"))
    tacsat$SI_DAY <- yday(tacsat$SI_DATIM)
  }
  if("SI_YEAR" %in% unique(c(by))){
    eflalo$SI_YEAR<- year(as.POSIXct(eflalo$LE_CDAT,format="%d/%m/%Y"))
    tacsat$SI_YEAR<- year(tacsat$SI_DATIM)
  }
  if("SI_MONTH"%in% unique(c(by))){
    eflalo$SI_MONTH<- month(as.POSIXct(eflalo$LE_CDAT,format="%d/%m/%Y"))
    tacsat$SI_MONTH<- month(tacsat$SI_DATIM)
  }
  if("SI_WEEK"%in% unique(c(by))){
    eflalo$SI_WEEK<- week(as.POSIXct(eflalo$LE_CDAT,format="%d/%m/%Y"))
    tacsat$SI_WEEK<- week(tacsat$SI_DATIM)
  }
  if("SI_QUARTER"%in% unique(c(by))){
    eflalo$SI_QUARTER<- quarter(as.POSIXct(eflalo$LE_CDAT,format="%d/%m/%Y"))
    tacsat$SI_QUARTER<- quarter(tacsat$SI_DATIM)
  }

  #- Select time levels + order these
  timePos<- c("SI_DAY","SI_WEEK","SI_MONTH","SI_QUARTER","SI_YEAR")
  byTime <- by[which(by %in% timePos)]
  byTime <- timePos[which(timePos %in% byTime)]

  #-----------------------------------------------------------------------------
  #-Add spatial location
  #-----------------------------------------------------------------------------
  if("LE_RECT" %in% by)
    tacsat$LE_RECT <- ICESrectangle(tacsat)
  if("LE_ICESAREA" %in% by){
    data(ICESareas)
    tacsat$LE_AREA <- ICESarea(tacsat,ICESareas)
    tacsat$LE_AREA[which(is.na(tacsat$LE_AREA)==T)] <- "OTHER"
    eflonlat       <- ICESrectangle2LonLat(eflalo$LE_RECT)
    eflalo$LE_AREA <- ICESarea(eflonlat,ICESareas)
    eflalo$LE_AREA[which(is.na(eflalo$LE_AREA)==T)] <- "OTHER"
  }

  #- Select area levels + order these
  areaPos<- c("LE_RECT","LE_ICESAREA")
  byArea <- by[which(by %in% areaPos)]
  byArea <- areaPos[which(areaPos %in% byArea)]

  #-----------------------------------------------------------------------------
  #- Calculate time possible to spend per day fishing
  #-----------------------------------------------------------------------------
  eflalo$LE_CDATIM <- as.POSIXct(eflalo$LE_CDAT,format="%d/%m/%Y",tz="GMT")
  eflalo$INTV      <- c(difftime(eflalo$FT_LDATIM,eflalo$FT_DDATIM,units="mins"))
  eflalo$FT_DURDAY <- ifelse(c(difftime(as.Date(eflalo$FT_LDATIM),as.Date(eflalo$FT_DDATIM),units="hours") == 0),
                        c(difftime(eflalo$FT_LDATIM,eflalo$FT_DDATIM,units="mins")),
                        ifelse(c(difftime(as.Date(eflalo$FT_DDATIM),as.Date(eflalo$LE_CDATIM),units="hours")==0),
                          c(difftime(eflalo$LE_CDATIM+(60*60*24),eflalo$FT_DDATIM,units="mins")),
                          ifelse(c(difftime(as.Date(eflalo$FT_LDATIM),as.Date(eflalo$LE_CDATIM),units="hours")==0),
                            c(difftime(eflalo$FT_LDATIM,eflalo$LE_CDATIM,units="mins")),
                            1440)))

  # Here there is still a problem because INTVDAY is calculated for catch days only, so you miss some effort of a whole trip
  eflalo$dummy     <- 1
  eflalo           <- merge(eflalo,aggregate(eflalo$dummy,by=list(eflalo$FT_REF,eflalo$LE_CDATIM),FUN=sum,na.rm=T),by.x=c("FT_REF","LE_CDATIM"),by.y=c("Group.1","Group.2"),all.x=T)
  colnames(eflalo)[length(colnames(eflalo))] <- "NR_FT_REF"
  if("SI_DAY" %in% by){
    eflalo$INTVDAY   <- eflalo$FT_DURDAY / eflalo$NR_FT_REF
  } else {
    eflalo$INTVDAY   <- eflalo$INTV / eflalo$NR_FT_REF
    }
  eflalo           <- eflalo[,-grep("dummy",colnames(eflalo))]
  eflalo           <- eflalo[,-grep("FT_DURDAY",colnames(eflalo))]
  eflalo           <- eflalo[,-grep("NR_FT_REF",colnames(eflalo))]

  #-----------------------------------------------------------------------------
  #- Check if all colums in both tacsat and eflalo are available
  #-----------------------------------------------------------------------------
  if(length(which(!by %in% colnames(tacsat) | !by %in% colnames(eflalo)))>0)
    stop("elements specified in 'by' are not available as columns in either tacsat or eflalo")

  #-----------------------------------------------------------------------------
  #- Estimate proportion fishing? (estimate on a basis which proportion of being out of harbour
  #  is used to fish and use that proportion to rescale eflalo effort)
  #-----------------------------------------------------------------------------

  subTacsat             <- subset(tacsat,SI_STATE != 0 & FT_REF != 0)
  subEflalo             <- eflalo

  if("SI_DAY" %in% by){
    subEflalo$INTVTRIP  <- subEflalo$INTV
    subEflalo$INTV      <- subEflalo$INTVDAY
  } else {
    subEflalo           <- subEflalo[!duplicated(subEflalo$FT_REF),]
  }
  if("LE_GEAR" %in% by & !"LE_GEAR" %in% colnames(subTacsat))
    subTacsat$LE_GEAR   <- subEflalo$LE_GEAR[match(subTacsat$FT_REF,subEflalo$FT_REF)]
  if(!"FT_REF" %in% colnames(subTacsat))
    stop("FT_REF column missing from tacsat")

  #- Calculate effort in tacsat and eflalo on a trip basis
  if(!"FT_REF" %in% by){
    by                  <- c(by,"FT_REF")
    warning("FT_REF has been added to 'by' statement")
  }
  byOrig                      <- by
  idxE                        <- 1:nrow(eflalo)
  idxEs                       <- 1:nrow(subEflalo)
  idxTs                       <- 1:nrow(subTacsat)

  #- First go: aggregate INTV on specified level and compare with each other
  if(length(by)>1){
    estEffTacsat              <- aggregate(subTacsat$INTV[idxTs],by=as.list(subTacsat[idxTs,by]),FUN=sum,na.rm=T)
    estEffEflalo              <- aggregate(subEflalo$INTV[idxEs],by=as.list(subEflalo[idxEs,by]),FUN=sum,na.rm=T)
  } else {
      estEffTacsat            <- aggregate(subTacsat$INTV[idxTs],by=list(subTacsat[idxTs,by]),FUN=sum,na.rm=T)
      estEffEflalo            <- aggregate(subEflalo$INTV[idxEs],by=list(subEflalo[idxEs,by]),FUN=sum,na.rm=T)
    }
  colnames(estEffTacsat)      <- c(by,"INTV")
  colnames(estEffEflalo)      <- c(by,"INTV")

  #- Reset to original estimPropFish
  mestEff                     <- merge(estEffTacsat,estEffEflalo,by=c(by))
  allGZero                    <- which(mestEff$INTV.x>0 & mestEff$INTV.y>0)
  by                          <- by[-grep("FT_REF",by)]

  #- Calculate estimated proportion of fishing for each eflalo ping
  if(length(by)>1){
    estimFish                 <- aggregate(log(mestEff$INTV.x[allGZero] / mestEff$INTV.y[allGZero]),
                                           by=as.list(mestEff[allGZero,by]),
                                           FUN=mean,na.rm=T)
  } else {
    estimFish                 <- aggregate(log(mestEff$INTV.x[allGZero] / mestEff$INTV.y[allGZero]),
                                           by=list(mestEff[allGZero,by]),
                                           FUN=mean,na.rm=T)
  }
  estimFish$x                 <- exp(estimFish$x)
  colnames(estimFish)         <- c(by,"x")
  allGZero                    <- which(mestEff$INTV.x>0 & mestEff$INTV.y>0)
  allMean                     <- exp(mean(log(mestEff$INTV.x[allGZero] / mestEff$INTV.y[allGZero])))

  #- Add fishing proportion parameter to eflalo (x)
  eflalo$x[idxE]              <- merge(eflalo[   idxE, ],estimFish,by=by,all.x=T)$x
  subEflalo$x[idxEs]          <- merge(subEflalo[idxEs,],estimFish,by=by,all.x=T)$x
  idxE                        <- which(is.na(   eflalo$x)==T)
  idxEs                       <- which(is.na(subEflalo$x)==T)
  idxTs                       <- which(!apply(t(t(subTacsat[,by])),1,paste,collapse="_") %in%
                                        apply(t(t(eflalo[which(is.na(eflalo$x)==F),by])),1,paste,collapse="_"))

  #- Second go: aggregate INTV on specified level and compare with each other but drop one element time at time
  if(length(idxTs)>0 & length(idxEs)>0){
    if(length(by)>1){
      for(i in 1:(length(by)-1)){
        by                      <- rev(rev(byOrig[-grep("FT_REF",byOrig)])[-c(1:i)])
        by                      <- c(by,"FT_REF")

        if(!"SI_DAY" %in% by){
          subEflalo             <- subEflalo[!duplicated(subEflalo$FT_REF),]
          subEflalo$INTV        <- subEflalo$INTVTRIP
        }

        if(length(by)>1){
          estEffTacsat          <- aggregate(subTacsat$INTV[idxTs],by=as.list(subTacsat[idxTs,by]),FUN=sum,na.rm=T)
          estEffEflalo          <- aggregate(subEflalo$INTV[idxEs],by=as.list(subEflalo[idxEs,by]),FUN=sum,na.rm=T)
        } else {
          estEffTacsat          <- aggregate(subTacsat$INTV[idxTs],by=list(subTacsat[idxTs,by]),FUN=sum,na.rm=T)
          estEffEflalo          <- aggregate(subEflalo$INTV[idxEs],by=list(subEflalo[idxEs,by]),FUN=sum,na.rm=T)
        }
        colnames(estEffEflalo)  <- c(by,"INTV")
        colnames(estEffTacsat)  <- c(by,"INTV")

        mestEff                 <- merge(estEffTacsat[,colnames(estEffEflalo)],estEffEflalo,by=c(by))
        allGZero                <- which(mestEff$INTV.x>0 & mestEff$INTV.y>0)

        #- Reset to original estimPropFish
        by                      <- by[-grep("FT_REF",by)]
        if(length(by)>1){
          estimFish             <- aggregate(log(mestEff$INTV.x[allGZero] / mestEff$INTV.y[allGZero]),
                                             by=as.list(mestEff[allGZero,by]),
                                             FUN=mean,na.rm=T)
        } else {
          estimFish             <- aggregate(log(mestEff$INTV.x[allGZero] / mestEff$INTV.y[allGZero]),
                                             by=list(mestEff[allGZero,by]),
                                             FUN=mean,na.rm=T)
        }
        estimFish$x             <- exp(estimFish$x)
        colnames(estimFish)     <- c(by,"x")
        eflalo$x[idxE]          <- merge(eflalo   [idxE, -grep("x",colnames(   eflalo))],estimFish,by=by,all.x=T)$x
        subEflalo$x[idxEs]      <- merge(subEflalo[idxEs,-grep("x",colnames(subEflalo))],estimFish,by=by,all.x=T)$x
        idxE                    <- which(is.na(   eflalo$x)==T)
        idxEs                   <- which(is.na(subEflalo$x)==T)
        idxTs                   <- which(!apply(t(t(subTacsat[,by])),1,paste,collapse="_") %in%
                                          apply(t(t(eflalo[which(is.na(eflalo$x)==F),by])),1,paste,collapse="_"))
      }
    }
  }
  #- Third go: aggregate INTV on FT_REF level only
  if(length(idxEs)>0 & length(idxTs)>0){
    if(!"SI_DAY" %in% by){
      subEflalo               <- subEflalo[!duplicated(subEflalo$FT_REF),]
      subEflalo$INTV          <- subEflalo$INTVTRIP
    }
    estEffTacsat              <- aggregate(subTacsat$INTV[idxTs],by=list(subTacsat[idxTs,"FT_REF"]),FUN=sum,na.rm=T)
    estEffEflalo              <- aggregate(subEflalo$INTV[idxEs],by=list(subEflalo[idxEs,"FT_REF"]),FUN=sum,na.rm=T)
    colnames(estEffEflalo)    <- c("FT_REF","INTV")
    colnames(estEffTacsat)    <- c("FT_REF","INTV")
    mestEff                   <- merge(estEffTacsat[,colnames(estEffEflalo)],estEffEflalo,by=c("FT_REF"))
    allGZero                  <- which(mestEff$INTV.x>0 & mestEff$INTV.y>0)
    estimFish                 <- aggregate(log(mestEff$INTV.x[allGZero] / mestEff$INTV.y[allGZero]),
                                  by=list(mestEff[allGZero,"FT_REF"]),
                                  FUN=mean,na.rm=T)
    estimFish$x               <- exp(estimFish$x)
    colnames(estimFish)       <- c("FT_REF","x")
    eflalo$x[idxE]            <- merge(eflalo   [idxE, -grep("x",colnames(   eflalo))],estimFish,by="FT_REF",all.x=T)$x
    subEflalo$x[idxEs]        <- merge(subEflalo[idxEs,-grep("x",colnames(subEflalo))],estimFish,by="FT_REF",all.x=T)$x
    idxE                      <- which(is.na(   eflalo$x)==T)
    idxEs                     <- which(is.na(subEflalo$x)==T)
    idxTs                     <- which(!apply(t(t(subTacsat[,by])),1,paste,collapse="_") %in%
                                        apply(t(t(eflalo[which(is.na(eflalo$x)==F),by])),1,paste,collapse="_"))
  }
  if(length(idxEs)>0){
    subEflalo$x[idxEs]        <- allMean
    eflalo$x[idxE]            <- allMean
    colnames(subEflalo)[length(colnames(subEflalo))] <- "PropFish"
    colnames(eflalo)[   length(colnames(   eflalo))] <- "PropFish"
  }
  colnames(eflalo)[   length(colnames(   eflalo))] <- "PropFish"
return(eflalo)}
