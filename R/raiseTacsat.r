raiseTacsat <- function(tacsat,eflalo,by=c("LE_GEAR","VE_REF","SI_DAY","LE_RECT","SI_YEAR"),sortBy=TRUE){


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
  
  #- If estimated proportion of fishing column is missing
  if(!"PropFish" %in% colnames(eflalo)){
    eflalo$PropFish             <- 1
  }

  #- Start to aggregate and raise tacsat, but drop elements through a for loop
  #  First aggregation without dropping any element
  tacsat$INTVR  <- 0

  cat("Start minutes eflalo",sum(eflalo$INTVDAY * eflalo$PropFish,na.rm=T),"\n")
  cat("Start minutes tacsat",sum(tacsat$INTV,na.rm=T),"\n")

  #- Order elements by uniqueness
  if(sortBy){
    if(length(by)>1)
      by          <- names(sort(apply(eflalo[,by],2,function(x){length(unique(x))})))
  }

  #- Drop element by element and merge
  byOrig    <- by
  INTVR     <- numeric(nrow(tacsat))
  INTVDAY   <- eflalo$INTVDAY
  for(j in 0:(length(by)-1)){
    if(j != 0)
      by  <- rev(rev(byOrig)[-c(1:j)])
    unEf  <- apply(unique(t(t(eflalo[,by]))),1,paste,collapse="_")
    pEf   <- apply(t(t(eflalo[,by])),1,paste,collapse="_")
    pTa   <- apply(t(t(tacsat[,by])),1,paste,collapse="_")
    cat("Raising by",by,"\n")
    for(i in 1:length(unEf)){
      idxT <- which(pTa %in% unEf[i])
      idxE <- which(pEf %in% unEf[i])
      if(length(idxT)>0 & length(idxE)>0){
        INTVR[idxT]   <- sum(INTVDAY[idxE] * eflalo$PropFish[idxE],na.rm=T)/sum(tacsat$INTV[idxT],na.rm=T) * tacsat$INTV[idxT] + INTVR[idxT]
        INTVDAY[idxE] <- 0
      }
    }
  }
  tacsat$INTVR        <- INTVR
  eflalo$INTVDAY      <- INTVDAY

  cat("Final minutes tacsat",sum(tacsat$INTVR,na.rm=T),"\n")
  cat("final remaing minutes eflalo",sum(eflalo$INTVDAY,na.rm=T),"\n")

return(tacsat)}

