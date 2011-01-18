splitAmongPings <- function(tacsat,eflalo,variable="all",level="day",conserve=T){

  #level: day,ICESrectangle,trip
  #conserve: T,F
  #variable: kgs,value,effort

  #- Create extra columns with time stamps
if(!"FT_REF" %in% colnames(tacsat)) stop("tacsat file needs FT_REF detailing trip number")
if(!"SI_STATE" %in% colnames(tacsat)) stop("tacsat file needs SI_STATE detailing activity of vessel")
if(level == "trip" & conserve == T) stop("conserve catches only at level = ICESrectangle or day")

if(!"SI_DATIM" %in%   colnames(tacsat)) tacsat$SI_DATIM     <- as.POSIXct(paste(tacsat$SI_DATE,   tacsat$SI_TIME,     sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
if(!"LE_CDATIM" %in%  colnames(eflalo)) eflalo$LE_CDATIM    <- as.POSIXct(eflalo$LE_CDAT,                                     tz="GMT", format="%d/%m/%Y")


  #- Levels have hierachy, and need to be suplemented with lower levels
if(level == "day"){               level <- c("day","ICESrectangle","trip")
} else {
    if(level == "ICESrectangle"){ level <- c("ICESrectangle","trip")
    } else {
        if(level == "trip"){
                                  level <- c("trip")
        }
      }
  }

  #- Add ID to keep track of merged and non-merged sets
tacsat$ID               <- 1:nrow(tacsat)

  #- identifyers of eflalo colnames
eflaloCol               <- colnames(eflalo)
kgs                     <- grep("KG",colnames(eflalo))
eur                     <- grep("EURO",colnames(eflalo))
kgeur                   <- function(x){return(c(grep("KG",x),grep("EURO",x)))}

  #- Subset tacsat file
tacsat                  <- subset(tacsat,SI_STATE != 0) #only attribute variable to fishing pings
tacsatTrip              <- subset(tacsat,FT_REF != 0)
remainTacsat            <- sort(unique(tacsatTrip$ID))

  #- Subset eflalo file
eflaloTrip              <- subset(eflalo,       FT_REF %in% sort(unique(tacsatTrip$FT_REF)))
eflaloNoTrip            <- subset(eflalo,      !FT_REF %in% sort(unique(tacsatTrip$FT_REF)))
eflaloVessel            <- subset(eflaloNoTrip, VE_REF %in% sort(unique(tacsatTrip$VE_REF)))
eflaloNoVessel          <- subset(eflaloNoTrip,!VE_REF %in% sort(unique(tacsatTrip$VE_REF)))

#-------------------------------------------------------------------------------
# 1a) Merge eflalo to tacsat with matching FT_REF
#-------------------------------------------------------------------------------

if("day" %in% level){
  if(!"SI_YEAR" %in% colnames(tacsatTrip))  tacsatTrip$SI_YEAR    <- an(format(tacsatTrip$SI_DATIM,format="%Y"))
  if(!"SI_DAY" %in%  colnames(tacsatTrip))  tacsatTrip$SI_DAY     <- an(format(tacsatTrip$SI_DATIM,format="%j"))
  if(!"LE_RECT" %in% colnames(tacsatTrip))  tacsatTrip$LE_RECT    <- ICESrectangle(tacsatTrip)

  if(!"SI_YEAR" %in% colnames(eflaloTrip))  eflaloTrip$SI_YEAR    <- an(format(eflaloTrip$LE_CDATIM,format="%Y"))
  if(!"SI_DAY" %in%  colnames(eflaloTrip))  eflaloTrip$SI_DAY     <- an(format(eflaloTrip$LE_CDATIM,format="%j"))

    #- Count pings in tacsat set
  nPings                <- countPings(~year+VE_REF+FT_REF+icesrectangle+day,tacsatTrip)

    #- Do the merging of eflalo to tacsat
  res           <- eflalo2Pings(eflaloTrip,tacsatTrip,nPings,c("SI_YEAR","VE_REF","FT_REF","LE_RECT","SI_DAY"),eflaloCol[c(kgs,eur)],remainTacsat)
  eflaloTrip    <- res[["eflalo"]]
  byDayTacsat   <- res[["tacsat"]]
  remainTacsat  <- res[["remainTacsat"]]
}
if("ICESrectangle" %in% level){
  if(!"SI_YEAR" %in% colnames(tacsatTrip))  tacsatTrip$SI_YEAR    <- an(format(tacsatTrip$SI_DATIM,format="%Y"))
  if(!"LE_RECT" %in% colnames(tacsatTrip))  tacsatTrip$LE_RECT    <- ICESrectangle(tacsatTrip)

  if(!"SI_YEAR" %in% colnames(eflaloTrip))  eflaloTrip$SI_YEAR    <- an(format(eflaloTrip$LE_CDATIM,format="%Y"))

    #- Count pings in tacsat set
  nPings                <- countPings(~year+VE_REF+FT_REF+icesrectangle,tacsatTrip)

    #- Do the merging of eflalo to tacsat
  res           <- eflalo2Pings(eflaloTrip,tacsatTrip,nPings,c("SI_YEAR","VE_REF","FT_REF","LE_RECT"),        eflaloCol[c(kgs,eur)],remainTacsat)
  eflaloTrip    <- res[["eflalo"]]
  byRectTacsat  <- res[["tacsat"]]
  remainTacsat  <- res[["remainTacsat"]]
}
if("trip" %in% level){
  if(!"SI_YEAR" %in% colnames(tacsatTrip))  tacsatTrip$SI_YEAR    <- an(format(tacsatTrip$SI_DATIM,format="%Y"))
  if(!"SI_YEAR" %in% colnames(eflaloTrip))  eflaloTrip$SI_YEAR    <- an(format(eflaloTrip$LE_CDATIM,format="%Y"))

    #- Count pings in tacsat set
  nPings                <- countPings(~year+VE_REF+FT_REF,tacsatTrip)
  
    #- Do the merging of eflalo to tacsat
  res           <- eflalo2Pings(eflaloTrip,tacsatTrip,nPings,c("SI_YEAR","VE_REF","FT_REF"),                  eflaloCol[c(kgs,eur)],remainTacsat)
  eflaloTrip    <- res[["eflalo"]]
  byTripTacsat  <- res[["tacsat"]]
  remainTacsat  <- res[["remainTacsat"]]
}
#-------------------------------------------------------------------------------
# 1b) Bind all tacsat files with matching FT_REF
#-------------------------------------------------------------------------------

if(length(remainTacsat) > 0) warning("Not all tacsat records with tripnumber have been merged!!")

if("day" %in% level){ tacsatFTREF <- rbind(byDayTacsat,byRectTacsat,byTripTacsat)
} else {
    if("ICESrectangle" %in% level){ tacsatFTREF <- rbind(byRectTacsat,byTripTacsat)
    } else { tacsatFTREF <- byTripTacsat}}
tacsatFTREF[,kgeur(colnames(tacsatFTREF))]    <- sweep(tacsatFTREF[,kgeur(colnames(tacsatFTREF))],1,tacsatFTREF$pings,"/")
tacsatFTREF$ID                                <- af(ac(tacsatFTREF$ID.x))
DT                                            <- data.table(tacsatFTREF)
eq1                                           <- c.listquote(paste("sum(",colnames(tacsatFTREF[,kgeur(colnames(tacsatFTREF))]),")",sep=""))
tacsatFTREF                                   <- DT[,eval(eq1),by=ID.x]; tacsatFTREF <- data.frame(tacsatFTREF); colnames(tacsatFTREF) <- c("ID",colnames(eflaloTrip[,kgeur(colnames(eflaloTrip))]))


#-------------------------------------------------------------------------------
# 2a) Merge eflalo to tacsat with no matching FT_REF
#-------------------------------------------------------------------------------

  #- If you don't want to loose catch or value data, conserve the non-merged
  #   eflalo catches and distribute these over the tacsat records
if(conserve == T){

  #-------------------------------------------------------------------------------
  # 2a-1) Merge eflalo to tacsat with matching VE_REF
  #-------------------------------------------------------------------------------

  if("day" %in% level){
    if(!"SI_YEAR" %in% colnames(tacsatTrip))  tacsatTrip$SI_YEAR    <- an(format(tacsatTrip$SI_DATIM,format="%Y"))
    if(!"SI_DAY" %in%  colnames(tacsatTrip))  tacsatTrip$SI_DAY     <- an(format(tacsatTrip$SI_DATIM,format="%j"))
    if(!"LE_RECT" %in% colnames(tacsatTrip))  tacsatTrip$LE_RECT    <- ICESrectangle(tacsatTrip)

    if(!"SI_YEAR" %in% colnames(eflaloVessel))  eflaloVessel$SI_YEAR    <- an(format(eflaloVessel$LE_CDATIM,format="%Y"))
    if(!"SI_DAY" %in%  colnames(eflaloVessel))  eflaloVessel$SI_DAY     <- an(format(eflaloVessel$LE_CDATIM,format="%j"))

      #- Count pings in tacsat set
    nPings                <- countPings(~year+VE_REF+icesrectangle+day,tacsatTrip)

      #- Do the merging of eflalo to tacsat
    res           <- eflalo2Pings(eflaloVessel,tacsatTrip,nPings,c("SI_YEAR","VE_REF","LE_RECT","SI_DAY"),      eflaloCol[c(kgs,eur)],NULL)
    eflaloVessel  <- res[["eflalo"]]
    byDayTacsat   <- res[["tacsat"]]
  }
  
  if("ICESrectangle" %in% level){
    if(!"SI_YEAR" %in% colnames(tacsatTrip))  tacsatTrip$SI_YEAR    <- an(format(tacsatTrip$SI_DATIM,format="%Y"))
    if(!"LE_RECT" %in% colnames(tacsatTrip))  tacsatTrip$LE_RECT    <- ICESrectangle(tacsatTrip)
    
    if(!"SI_YEAR" %in% colnames(eflaloVessel))  eflaloVessel$SI_YEAR    <- an(format(eflaloVessel$LE_CDATIM,format="%Y"))

      #- Count pings in tacsat set
    nPings                <- countPings(~year+VE_REF+icesrectangle,tacsatTrip)

      #- Do the merging of eflalo to tacsat
    res           <- eflalo2Pings(eflaloVessel,tacsatTrip,nPings,c("SI_YEAR","VE_REF","LE_RECT"),               eflaloCol[c(kgs,eur)],NULL)
    eflaloVessel  <- res[["eflalo"]]
    byRectTacsat  <- res[["tacsat"]]
  }
  if(TRUE){ #-For remainder of vessel merging not at ICESrectangle level
    if(!"SI_YEAR" %in% colnames(tacsatTrip))  tacsatTrip$SI_YEAR    <- an(format(tacsatTrip$SI_DATIM,format="%Y"))
    if(!"SI_YEAR" %in% colnames(eflaloVessel))  eflaloVessel$SI_YEAR    <- an(format(eflaloVessel$LE_CDATIM,format="%Y"))

      #- Count pings in tacsat set
    nPings                <- countPings(~year+VE_REF,tacsatTrip)

      #- Do the merging of eflalo to tacsat
    res           <- eflalo2Pings(eflaloVessel,tacsatTrip,nPings,c("SI_YEAR","VE_REF" ),               eflaloCol[c(kgs,eur)],NULL)
    eflaloVessel  <- res[["eflalo"]]
    byVessTacsat  <- res[["tacsat"]]
  }
  
  #-------------------------------------------------------------------------------
  # 2b-1) Bind all tacsat files with matching VE_REF
  #-------------------------------------------------------------------------------

  if("day" %in% level){ tacsatVEREF <- rbind(byDayTacsat,byRectTacsat,byVessTacsat)
  } else {
      if("ICESrectangle" %in% level){
        tacsatVEREF <- rbind(byRectTacsat,byVessTacsat)
      } else { tacsatVEREF <- byVessTacsat}}
  tacsatVEREF[,kgeur(colnames(tacsatVEREF))]    <- sweep(tacsatVEREF[,kgeur(colnames(tacsatVEREF))],1,tacsatVEREF$pings,"/")
  tacsatVEREF$ID                                <- af(ac(tacsatVEREF$ID.x))
  DT                                            <- data.table(tacsatVEREF)
  eq1                                           <- c.listquote(paste("sum(",colnames(tacsatVEREF[,kgeur(colnames(tacsatVEREF))]),")",sep=""))
  tacsatVEREF                                   <- DT[,eval(eq1),by=ID.x]; tacsatVEREF <- data.frame(tacsatVEREF); colnames(tacsatVEREF) <- c("ID",colnames(eflaloVessel[,kgeur(colnames(eflaloVessel))]))

  #-------------------------------------------------------------------------------
  # 2a-2) Merge eflalo to tacsat with no matching FT_REF or VE_REF
  #-------------------------------------------------------------------------------
  if("day" %in% level){
    if(!"SI_YEAR" %in% colnames(tacsatTrip))  tacsatTrip$SI_YEAR    <- an(format(tacsatTrip$SI_DATIM,format="%Y"))
    if(!"SI_DAY" %in%  colnames(tacsatTrip))  tacsatTrip$SI_DAY     <- an(format(tacsatTrip$SI_DATIM,format="%j"))
    if(!"LE_RECT" %in% colnames(tacsatTrip))  tacsatTrip$LE_RECT    <- ICESrectangle(tacsatTrip)

    if(!"SI_YEAR" %in% colnames(eflaloNoVessel))  eflaloNoVessel$SI_YEAR    <- an(format(eflaloNoVessel$LE_CDATIM,format="%Y"))
    if(!"SI_DAY" %in%  colnames(eflaloNoVessel))  eflaloNoVessel$SI_DAY     <- an(format(eflaloNoVessel$LE_CDATIM,format="%j"))

      #- Count pings in tacsat set
    nPings                <- countPings(~year+icesrectangle+day,tacsatTrip)

     #- Do the merging of eflalo to tacsat
    res               <- eflalo2Pings(eflaloNoVessel,tacsatTrip,nPings,c("SI_YEAR","LE_RECT","SI_DAY"),               eflaloCol[c(kgs,eur)],NULL)
    eflaloNoVessel    <- res[["eflalo"]]
    byDayTacsat       <- res[["tacsat"]]
  }

  if("ICESrectangle" %in% level){
    if(!"SI_YEAR" %in% colnames(tacsatTrip))  tacsatTrip$SI_YEAR    <- an(format(tacsatTrip$SI_DATIM,format="%Y"))
    if(!"LE_RECT" %in% colnames(tacsatTrip))  tacsatTrip$LE_RECT    <- ICESrectangle(tacsatTrip)

    if(!"SI_YEAR" %in% colnames(eflaloNoVessel))  eflaloNoVessel$SI_YEAR    <- an(format(eflaloNoVessel$LE_CDATIM,format="%Y"))

      #- Count pings in tacsat set
    nPings            <- countPings(~year+icesrectangle,tacsatTrip)

     #- Do the merging of eflalo to tacsat
    res               <- eflalo2Pings(eflaloNoVessel,tacsatTrip,nPings,c("SI_YEAR","LE_RECT"),                        eflaloCol[c(kgs,eur)],NULL)
    eflaloNoVessel    <- res[["eflalo"]]
    byRectTacsat      <- res[["tacsat"]]
  }
  if(TRUE){ #-For remainder of merging not at ICESrectangle level
    if(!"SI_YEAR" %in% colnames(tacsatTrip))  tacsatTrip$SI_YEAR    <- an(format(tacsatTrip$SI_DATIM,format="%Y"))
    if(!"SI_YEAR" %in% colnames(eflaloNoVessel))  eflaloNoVessel$SI_YEAR    <- an(format(eflaloNoVessel$LE_CDATIM,format="%Y"))
    
      #- Count pings in tacsat set
    nPings            <- countPings(~year,tacsatTrip)

     #- Do the merging of eflalo to tacsat
    res               <- eflalo2Pings(eflaloNoVessel,tacsatTrip,nPings,c("SI_YEAR"),                        eflaloCol[c(kgs,eur)],NULL)
    eflaloNoVessel    <- res[["eflalo"]]
    byVessTacsat      <- res[["tacsat"]]
  }
  #-------------------------------------------------------------------------------
  # 2b-2) Bind all tacsat files with no matching FT_REF or VE_REF
  #-------------------------------------------------------------------------------

  if("day" %in% level){ tacsatREF <- rbind(byDayTacsat,byRectTacsat,byVessTacsat)
  } else {
      if("ICESrectangle" %in% level){ tacsatREF <- rbind(byRectTacsat,byVessTacsat)
      } else { tacsatREF <- byVessTacsat}}
  tacsatREF[,kgeur(colnames(tacsatREF))]      <- sweep(tacsatREF[,kgeur(colnames(tacsatREF))],1,tacsatREF$pings,"/")
  tacsatREF$ID                                <- af(ac(tacsatREF$ID.x))
  DT                                          <- data.table(tacsatREF)
  eq1                                         <- c.listquote(paste("sum(",colnames(tacsatREF[,kgeur(colnames(tacsatREF))]),")",sep=""))
  tacsatREF                                   <- DT[,eval(eq1),by=ID.x]; tacsatREF <- data.frame(tacsatREF); colnames(tacsatREF) <- c("ID",colnames(eflaloVessel[,kgeur(colnames(eflaloVessel))]))
}#End conserve

#-------------------------------------------------------------------------------
# 3) Merge all tacsat files together and return
#-------------------------------------------------------------------------------

if(conserve==T){
  tacsatTot       <- rbind(tacsatFTREF,tacsatVEREF,tacsatREF)
  DT              <- data.table(tacsatTot)
  eq1             <- c.listquote(paste("sum(",colnames(tacsatTot[,kgeur(colnames(tacsatTot))]),")",sep=""))
  tacsatTot       <- DT[,eval(eq1),by=ID]; tacsatTot <- data.frame(tacsatTot); colnames(tacsatTot) <- c("ID",colnames(eflaloVessel[,kgeur(colnames(eflaloVessel))]))
  tacsatReturn    <- merge(tacsat,tacsatTot,by="ID",all.x=T)
  if(variable == "value") tacsatReturn <- tacsatReturn[,c(1:dim(tacsat)[2],grep("EURO",colnames(tacsatReturn)))]
  if(variable == "kgs")   tacsatReturn <- tacsatReturn[,c(1:dim(tacsat)[2],grep("KG",colnames(tacsatReturn)))]
  if(variable == "all")   tacsatReturn <- tacsatReturn
} else {
    tacsatReturn  <- merge(tacsat,tacsatFTREF,by="ID",all.x=T)
    if(variable == "value") tacsatReturn <- tacsatReturn[,c(1:dim(tacsat)[2],grep("EURO",colnames(tacsatReturn)))]
    if(variable == "kgs")   tacsatReturn <- tacsatReturn[,c(1:dim(tacsat)[2],grep("KG",colnames(tacsatReturn)))]
    if(variable == "all")   tacsatReturn <- tacsatReturn
  }

return(orderBy(~ID,data=tacsatReturn)[,-grep("ID",colnames(tacsatReturn))])}



