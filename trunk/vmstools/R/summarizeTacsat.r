summarizeTacsat <- function(tacsat){

  nrVessels     <- length(unique(tacsat$VE_REF))
  nrCountries   <- length(unique(tacsat$VE_COU))

  if (!"SI_DATIM" %in% colnames(tacsat))
        tacsat$SI_DATIM <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,
            sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M")

  if("FT_REF" %in% colnames(tacsat)){
    totalEffort  <- sum(intervalTacsat(tacsat,level="trip",fill.na=T)$INTV,na.rm=T)/60
  } else {
    totalEffort  <- sum(intervalTacsat(tacsat,level="vessel",fill.na=T)$INTV,na.rm=T)/60
  }
  return(as.data.frame(
          cbind(desc= c("nrCountries","nrVessels","minLon","maxLon","minLat","maxLat","minTime","maxTime","minHeading","maxHeading","minSpeed","maxSpeed","effort(hr)"),
                value=c(nrCountries,nrVessels,round(range(tacsat$SI_LONG,na.rm=T),3),round(range(tacsat$SI_LATI,na.rm=T),3),
                        ac(range(tacsat$SI_DATIM,na.rm=T)[1]),ac(range(tacsat$SI_DATIM,na.rm=T)[2]),
                        range(tacsat$SI_HE,na.rm=T),range(tacsat$SI_SP,na.rm=T),round(totalEffort,1))),stringsAsFactors=F))}
  
summarizeEflalo <- function(eflalo){

  nrVessels     <- length(unique(eflalo$VE_REF))
  nrCountries   <- length(unique(eflalo$VE_COU))

  if (!"FT_DDATIM" %in% colnames(eflalo))
        eflalo$FT_DDATIM <- as.POSIXct(paste(eflalo$FT_DDAT, eflalo$FT_DTIME,
            sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M")
  if (!"FT_LDATIM" %in% colnames(eflalo))
        eflalo$FT_LDATIM <- as.POSIXct(paste(eflalo$FT_LDAT, eflalo$FT_LTIME,
            sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M")

  totalEffort   <- sum(difftime(eflalo$FT_LDATIM,eflalo$FT_DDATIM,units="hours"),na.rm=T)
  totalKW       <- sum(eflalo$VE_KW,na.rm=T)
  meanLength    <- mean(eflalo$VE_LEN,na.rm=T)
  gears         <- names(rev(table(eflalo$LE_GEAR)))[1:3]
  catchvals     <- colSums(eflalo[,kgeur(colnames(eflalo))],na.rm=T)
  catches       <- catchvals[grep("LE_KG_",names(catchvals))]
  values        <- catchvals[grep("LE_EURO_",names(catchvals))]
  
  return(data.frame(
          cbind(desc= c("nrCountries","nrVessels","minTime","maxTime","effort(hr)","totalKW","meanLength",
                        "gear1","gear2","gear3","catch1","catch2","catch3","value1","value2","value3"),
                value=c(nrCountries,nrVessels,ac(range(eflalo$FT_DDATIM)[1]),ac(range(eflalo$FT_LDATIM)[2]),
                        round(totalEffort,1),round(totalKW,1),round(meanLength,1),gears,
                        paste(names(rev(sort(catches)))[1:3],round(rev(sort(catches))[1:3],0)),
                        paste(names(rev(sort(values )))[1:3],round(rev(sort(values ))[1:3],0)))),
                stringsAsFactors=F))}
  








