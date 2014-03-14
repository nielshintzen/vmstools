#-------------------------------------------------------------------------------
#
# Benthis WP2 workflow
#
# Designed by: Francois Bastardie, Niels Hintzen
# Runs with: R3.0.2 and R2.15.x
#
# VMStools version: 0.70
#
#-------------------------------------------------------------------------------

rm(list=ls())
library(vmstools)
library(maps)
library(mapdata)

if(.Platform$OS.type == "unix") {}
 codePath  <- file.path("~","BENTHIS")
 dataPath  <- file.path("~","BENTHIS","EflaloAndTacsat")
 outPath   <- file.path("~","BENTHIS", "outputs")
 polPath   <- file.path("~","BENTHIS", "BalanceMaps")
 
 
if(.Platform$OS.type == "windows") {
 codePath  <- "C:/BENTHIS/"
 dataPath  <- "C:/BENTHIS/EflaloAndTacsat/"
 outPath   <- "C:/BENTHIS/outputs/"
 polPath   <- "C:/BENTHIS/BalanceMaps"
 }


#a_year      <- 2010
#a_year      <- 2011
#a_year      <- 2012
a_year      <- 2013
dir.create(file.path(outPath))
dir.create(file.path(outPath, a_year))

if(TRUE){

  load(file.path(dataPath,paste("tacsat_", a_year,".RData", sep=''))); # get the tacsat object
  load(file.path(dataPath,paste("eflalo_", a_year,".RData", sep=''))); # get the eflalo object
  tacsat <- formatTacsat(tacsat) # format each of the columns to the specified class
  eflalo <- formatEflalo(eflalo) # format each of the columns to the specified class

  # drop the species catch and euro composition
  idx                   <- kgeur(colnames(eflalo))
  eflalo$LE_KG_SPECS    <- rowSums(eflalo[,grep("LE_KG_",colnames(eflalo))],na.rm=T)
  eflalo$LE_EURO_SPECS  <- rowSums(eflalo[,grep("LE_EURO_",colnames(eflalo))],na.rm=T)
  eflalo                <- eflalo[,-idx]

  # country-specific
  ctry   <- "DNK"
  eflalo <- eflalo[ grep(ctry, as.character(eflalo$VE_REF)),]  # keep the national vessels only.
  #VMS_ping_rate_in_hour <- 115/60 # Dutch data (rev(sort(table(intervalTacsat(sortTacsat(tacsat),level="vessel")$INTV))[1])
  VMS_ping_rate_in_hour <- 1 # e.g. 1 hour for Denmark (rev(sort(table(intervalTacsat(sortTacsat(tacsat),level="vessel")$INTV))[1])
  
  # Gear codes to keep (with assumed severe bottom impact)
  gears2keep            <- c("TBB","OTT","OTB","SSC","SDN","PTB","DRB")
  towedGears            <- c("TBB","OTT","OTB","PTB","DRB")
  seineGears            <- c("SSC","SDN")

  if(.Platform$OS.type == "windows")
    data(euharbours)
  if(.Platform$OS.type == "unix")
    data(harbours)
  data(ICESareas)
  data(europa)

  #-----------------------------------------------------------------------------
  # A FUNCTION (FOR LATER USE)
  #-----------------------------------------------------------------------------
  compute_swept_area <- function(
                              tacsatIntGearVEREF=tacsatIntGearVEREF,
                              gear_param_per_metier=gear_param_per_metier,
                              towedGears=towedGears,
                              seineGears=seineGears,
                              VMS_ping_rate_in_hour=VMS_ping_rate_in_hour,
                              already_informed_width_for=NULL
                              ){

  if(is.null(already_informed_width_for)){
     tacsatIntGearVEREF <- tacsatIntGearVEREF[,!colnames(tacsatIntGearVEREF) %in%
                          c('GEAR_WIDTH', 'GEAR_WIDTH_LOWER', 'GEAR_WIDTH_UPPER', 'SWEPT_AREA_KM2', 'SWEPT_AREA_KM2_LOWER', 'SWEPT_AREA_KM2_UPPER')]  # remove columns if exists
     } else{
     tacsatIntGearVEREF <- tacsatIntGearVEREF[,!colnames(tacsatIntGearVEREF) %in%
                          c('SWEPT_AREA_KM2', 'SWEPT_AREA_KM2_LOWER', 'SWEPT_AREA_KM2_UPPER')]  # remove columns if exists
     
     }
  
  if(is.null(already_informed_width_for)){
  # MERGE WITH GEAR WIDTH
  GearWidth                   <- tacsatIntGearVEREF[!duplicated(data.frame(tacsatIntGearVEREF$VE_REF,tacsatIntGearVEREF$LE_MET,tacsatIntGearVEREF$VE_KW,tacsatIntGearVEREF$VE_LEN)), ]
  GearWidth                   <- GearWidth[,c('VE_REF','LE_MET','VE_KW', 'VE_LEN') ]
  GearWidth$GEAR_WIDTH        <- NA
  GearWidth$GEAR_WIDTH_LOWER  <- NA
  GearWidth$GEAR_WIDTH_UPPER  <- NA
  for (i in 1:nrow(GearWidth)) { # brute force...
    kW      <- GearWidth$VE_KW[i]
    LOA     <- GearWidth$VE_LEN[i]
    this    <- gear_param_per_metier[gear_param_per_metier$a_metier==as.character(GearWidth$LE_MET[i]),]
    a <- NULL ; b <- NULL
    a       <- this[this$param=='a', 'Estimate']
    b       <- this[this$param=='b', 'Estimate']
    GearWidth[i,"GEAR_WIDTH"]  <-   eval(parse(text= as.character(this[1, 'equ']))) / 1000 # converted in km
    a       <- this[this$param=='a', 'Estimate']
    b       <- this[this$param=='b', 'Estimate'] +2*this[this$param=='b', 'Std..Error']
    GearWidth[i,"GEAR_WIDTH_UPPER"]  <-  eval(parse(text= as.character(this[1, 'equ']))) / 1000 # converted in km
    a       <- this[this$param=='a', 'Estimate']
    b       <- this[this$param=='b', 'Estimate'] -2*this[this$param=='b', 'Std..Error']
    GearWidth[i,"GEAR_WIDTH_LOWER"]  <-  eval(parse(text= as.character(this[1, 'equ']))) / 1000 # converted in km
  }
  tacsatIntGearVEREF                    <- merge(tacsatIntGearVEREF, GearWidth,by=c("VE_REF","LE_MET","VE_KW","VE_LEN"),
                                              all.x=T,all.y=F)

  }
                                              
  #  the swept area (note that could work oustide the loop area as well....)
  # for the trawlers...
  if(tacsatIntGearVEREF$LE_GEAR[1] %in% towedGears){
        tacsatIntGearVEREF$SWEPT_AREA_KM2 <- NA
        tacsatIntGearVEREF <- orderBy(~SI_DATIM,data=tacsatIntGearVEREF)
        a_dist             <- distance(c(tacsatIntGearVEREF$SI_LONG[-1],0),  c(tacsatIntGearVEREF$SI_LATI[-1],0),
                                         tacsatIntGearVEREF$SI_LONG, tacsatIntGearVEREF$SI_LATI)
        a_dist[length(a_dist)] <- rev(a_dist)[2]
        tacsatIntGearVEREF$SWEPT_AREA_KM2 <- a_dist * tacsatIntGearVEREF$GEAR_WIDTH
        tacsatIntGearVEREF$SWEPT_AREA_KM2_LOWER <- a_dist * tacsatIntGearVEREF$GEAR_WIDTH_LOWER
        tacsatIntGearVEREF$SWEPT_AREA_KM2_UPPER <- a_dist * tacsatIntGearVEREF$GEAR_WIDTH_UPPER
  }

  # for the seiners...
  if(tacsatIntGearVEREF$LE_GEAR[1]  %in% seineGears){
      tacsatIntGearVEREF$SWEPT_AREA_KM2         <- pi*(tacsatIntGearVEREF$GEAR_WIDTH/(2*pi))^2
      tacsatIntGearVEREF$SWEPT_AREA_KM2_LOWER   <- pi*(tacsatIntGearVEREF$GEAR_WIDTH_LOWER/(2*pi))^2
      tacsatIntGearVEREF$SWEPT_AREA_KM2_UPPER   <- pi*(tacsatIntGearVEREF$GEAR_WIDTH_UPPER/(2*pi))^2

      haul_duration                           <- 3 # assumption of a mean duration based from questionnaires to seiners
      tacsatIntGearVEREF$SWEPT_AREA_KM2         <- tacsatIntGearVEREF$SWEPT_AREA_KM2 * VMS_ping_rate_in_hour / haul_duration # correction to avoid counting the same circle are several time.
      tacsatIntGearVEREF$SWEPT_AREA_KM2_LOWER   <- tacsatIntGearVEREF$SWEPT_AREA_KM2_LOWER * VMS_ping_rate_in_hour / haul_duration # correction to avoid counting the same circle are several time.
      tacsatIntGearVEREF$SWEPT_AREA_KM2_UPPER   <- tacsatIntGearVEREF$SWEPT_AREA_KM2_UPPER * VMS_ping_rate_in_hour / haul_duration # correction to avoid counting the same circle are several time.
      idx                                     <- grep('SSC', as.character(tacsatIntGearVEREF$LE_GEAR))
      tacsatIntGearVEREF[idx, 'SWEPT_AREA_KM2'] <- tacsatIntGearVEREF[idx, 'SWEPT_AREA_KM2'] *1.5 # ad hoc correction to account for the SSC specificities
      tacsatIntGearVEREF[idx, 'SWEPT_AREA_KM2_LOWER'] <- tacsatIntGearVEREF[idx, 'SWEPT_AREA_KM2_LOWER'] *1.5 # ad hoc correction to account for the SSC specificities
      tacsatIntGearVEREF[idx, 'SWEPT_AREA_KM2_UPPER'] <- tacsatIntGearVEREF[idx, 'SWEPT_AREA_KM2_UPPER'] *1.5 # ad hoc correction to account for the SSC specificities
   }

   return(tacsatIntGearVEREF)
   }

  
  #-----------------------------------------------------------------------------
  # Cleaning tacsat (keep track of removed records)
  #-----------------------------------------------------------------------------
  remrecsTacsat      <- matrix(NA,nrow=6,ncol=2,dimnames= list(c("total","duplicates","notPossible",
                                                                 "pseudoDuplicates","harbour","land"),
                                                               c("rows","percentage")))
  remrecsTacsat["total",] <- c(nrow(tacsat),"100%")

  # Remove duplicate records
  tacsat$SI_DATIM     <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "),
                                   tz="GMT", format="%d/%m/%Y  %H:%M")
  uniqueTacsat        <- paste(tacsat$VE_REF,tacsat$SI_LATI,tacsat$SI_LONG,tacsat$SI_DATIM)
  tacsat              <- tacsat[!duplicated(uniqueTacsat),]
  remrecsTacsat["duplicates",] <- c(nrow(tacsat),100+round((nrow(tacsat) -
                                   an(remrecsTacsat["total",1]))/an(remrecsTacsat["total",1])*100,2))

  # Remove points that cannot be possible
  spThres             <- 20   #Maximum speed threshold in analyses in nm
  idx                 <- which(abs(tacsat$SI_LATI) > 90 | abs(tacsat$SI_LONG) > 180)
  idx                 <- unique(c(idx,which(tacsat$SI_HE < 0 | tacsat$SI_HE > 360)))
  idx                 <- unique(c(idx,which(tacsat$SI_SP > spThres)))
  if(length(idx)>0)
    tacsat            <- tacsat[-idx,]
  remrecsTacsat["notPossible",] <- c(nrow(tacsat),100+round((nrow(tacsat) -
                                     an(remrecsTacsat["total",1]))/an(remrecsTacsat["total",1])*100,2))

  # Remove points which are pseudo duplicates as they have an interval rate < x minutes
  intThres            <- 5    # Minimum difference in time interval in minutes to prevent pseudo duplicates
  tacsat              <- sortTacsat(tacsat)
  tacsatp             <- intervalTacsat(tacsat,level="vessel",fill.na=T)
  tacsat              <- tacsatp[which(tacsatp$INTV > intThres | is.na(tacsatp$INTV)==T),-grep("INTV",colnames(tacsatp))]
  remrecsTacsat["pseudoDuplicates",] <- c(nrow(tacsat),100+round((nrow(tacsat) -
                                          an(remrecsTacsat["total",1]))/an(remrecsTacsat["total",1])*100,2))

  # Remove points in harbour 
  idx             <- pointInHarbour(tacsat$SI_LONG,tacsat$SI_LATI,harbours)
  pih             <- tacsat[which(idx == 1),]
  save(pih,file=paste(outPath, a_year, "pointInHarbour.RData",sep=""))
  tacsat          <- tacsat[which(idx == 0),]
  remrecsTacsat["harbour",] <- c(nrow(tacsat),100+round((nrow(tacsat) -
                  an(remrecsTacsat["total",1]))/an(remrecsTacsat["total",1])*100,2))

  # Remove points on land
  pols                <- lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),
                                                function(x){data.frame(SI_LONG=subset(europa,SID==x)$X,
                                                                       SI_LATI=subset(europa,SID==x)$Y)}))
  idx                 <- pointOnLand(tacsat,pols)
  pol                 <- tacsat[which(idx == 1),]
  save(pol,file=file.path(outPath,a_year,"pointOnLand.RData"))
  tacsat              <- tacsat[which(idx == 0),]
  remrecsTacsat["land",] <- c(nrow(tacsat),100+round((nrow(tacsat) -
                an(remrecsTacsat["total",1]))/an(remrecsTacsat["total",1])*100,2))

  # Save the remrecsTacsat file
  save(remrecsTacsat,file=file.path(outPath,a_year,"remrecsTacsat.RData"))

  # Save the cleaned tacsat file
  save(tacsat,file=file.path(outPath,a_year,"cleanTacsat.RData"))

  #-----------------------------------------------------------------------------
  # Cleaning Eflalo
  #-----------------------------------------------------------------------------

  # Keep track of removed points
  remrecsEflalo     <- matrix(NA,nrow=5,ncol=2,dimnames=list(c("total","duplicated","impossible time",
                                                               "before 1st Jan","departArrival"),
                                                             c("rows","percentage")))
  remrecsEflalo["total",] <- c(nrow(eflalo),"100%")

  # Remove non-unique trip numbers
  eflalo            <- eflalo[!duplicated(paste(eflalo$LE_ID,eflalo$LE_CDAT,sep="-")),]
  remrecsEflalo["duplicated",] <- c(nrow(eflalo),100+round((nrow(eflalo) -
                                    an(remrecsEflalo["total",1]))/an(remrecsEflalo["total",1])*100,2))

  # Remove impossible time stamp records
  eflalo$FT_DDATIM  <- as.POSIXct(paste(eflalo$FT_DDAT,eflalo$FT_DTIME, sep = " "),
                               tz = "GMT", format = "%d/%m/%Y  %H:%M")
  eflalo$FT_LDATIM  <- as.POSIXct(paste(eflalo$FT_LDAT,eflalo$FT_LTIME, sep = " "),
                               tz = "GMT", format = "%d/%m/%Y  %H:%M")

  eflalo            <- eflalo[!(is.na(eflalo$FT_DDATIM) |is.na(eflalo$FT_LDATIM)),]
  remrecsEflalo["impossible time",] <- c(nrow(eflalo),100+round((nrow(eflalo) -
                  an(remrecsEflalo["total",1]))/an(remrecsEflalo["total",1])*100,2))

  # Remove trip starting before 1st Jan
  year              <- min(year(eflalo$FT_DDATIM))
  eflalo            <- eflalo[eflalo$FT_DDATIM>=strptime(paste(year,"-01-01 00:00:00",sep=''),
                                                             "%Y-%m-%d %H:%M"),]
  remrecsEflalo["before 1st Jan",] <- c(nrow(eflalo),100+round((nrow(eflalo) -
                  an(remrecsEflalo["total",1]))/an(remrecsEflalo["total",1])*100,2))

  # Remove records with arrival date before departure date
  eflalop           <- eflalo
  eflalop$FT_DDATIM <- as.POSIXct(paste(eflalo$FT_DDAT,  eflalo$FT_DTIME,   sep=" "),
                                tz="GMT", format="%d/%m/%Y  %H:%M")
  eflalop$FT_LDATIM <- as.POSIXct(paste(eflalo$FT_LDAT,  eflalo$FT_LTIME,   sep=" "),
                                tz="GMT", format="%d/%m/%Y  %H:%M")
  idx               <- which(eflalop$FT_LDATIM >= eflalop$FT_DDATIM)
  eflalo            <- eflalo[idx,]
  remrecsEflalo["departArrival",] <- c(nrow(eflalo),100+round((nrow(eflalo) -
                  an(remrecsEflalo["total",1]))/an(remrecsEflalo["total",1])*100,2))

  # Save the remrecsEflalo file
  save(remrecsEflalo,file=file.path(outPath,a_year,"remrecsEflalo.RData"))

  # Save the cleaned eflalo file
  save(eflalo,file=file.path(outPath,a_year,"cleanEflalo.RData"))

  #-----------------------------------------------------------------------------
  # Make gear code selection and calculate effort for each gear
  #-----------------------------------------------------------------------------

  # effort < 15m vs >15m
  eflalo$length_class   <- cut(as.numeric(as.character(eflalo$VE_LEN)), breaks=c(0,15,100))   #DCF but VMS!

  # compute effort
  eflalo                <- subset(eflalo,FT_REF != 0)
  eflalo                <- orderBy(~VE_REF+FT_DDATIM+FT_REF, data=eflalo)
  eflalo$ID             <- paste(eflalo$VE_REF,eflalo$FT_REF,sep="")
  eflalo$LE_EFF         <- an(difftime(eflalo$FT_LDATIM, eflalo$FT_DDATIM, units="hours"))
  eflalo$dummy          <- 1
  eflalo$LE_EFF         <- eflalo$LE_EFF / merge(eflalo,aggregate(eflalo$dummy,by=list(eflalo$ID),FUN=sum),by.x="ID",by.y="Group.1",all.x=T)$x
  eflalo                <- eflalo[which(eflalo$LE_GEAR %in% gears2keep),] 
  aggregate(eflalo$LE_EFF, list(eflalo$length_class), sum,na.rm=T)

  gc(reset=TRUE)

  #-----------------------------------------------------------------------------
  # Merge eflalo and tacsat
  #-----------------------------------------------------------------------------

  tacsatp               <- mergeEflalo2Tacsat(eflalo,tacsat)

  tacsatp$LE_GEAR       <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]
  tacsatp$VE_LEN        <- eflalo$VE_LEN[ match(tacsatp$FT_REF,eflalo$FT_REF)]
  tacsatp$LE_MET        <- eflalo$LE_MET[ match(tacsatp$FT_REF,eflalo$FT_REF)]
  tacsatp$VE_KW         <- eflalo$VE_KW[ match(tacsatp$FT_REF,eflalo$FT_REF)]
  if("LE_WIDTH" %in% colnames(eflalo))
    tacsatp$LE_WIDTH    <- eflalo$LE_WIDTH[ match(tacsatp$FT_REF,eflalo$FT_REF)]
  save(tacsatp,file=file.path(outPath,a_year,"tacsatMerged.RData"))

  # Save not merged tacsat data
  tacsatpmin            <- subset(tacsatp,FT_REF == 0)
  save(tacsatpmin, file=file.path(outPath,a_year,"tacsatNotMerged.RData"))

  tacsatp               <- subset(tacsatp,FT_REF != 0)
  
  
   #-----------------------------------------------------------------------------
  # transform into WP2 BENTHIS metier - ADAPT TO YOUR OWN METIER LIST!!!
  #-----------------------------------------------------------------------------
   ctry <- 'DNK'
   if(ctry=="DNK"){
   tacsatp$LE_MET_init    <- tacsatp$LE_MET
   tacsatp$LE_MET         <- factor(tacsatp$LE_MET)            
   print(levels(tacsatp$LE_MET))
   if(a_year=="2010"){
     levels(tacsatp$LE_MET) <-   c(  ## REPLACE LEVELS WITH CAUTION ## adapt to your own list!!
     "DRB_MOL", "NA",  "OT_CRU","OT_CRU", "OT_CRU",   "OT_CRU",  "OT_CRU",  
     "OT_CRU",  "OT_SPF",     "OT_DMF", "OT_DMF", "OT_DMF",   "OT_DMF", "OT_SPF",  
     "OT_MIX_DMF_PEL",   "OT_DMF",  "OT_DMF",   "OT_DMF",  "OT_DMF",  "OT_SPF",  "OT_SPF",  
     "OT_SPF",  "OT_SPF",   "OT_SPF",   "OT_SPF",   "OT_CRU",  "OT_SPF",     "OT_DMF",
     "OT_DMF",   "OT_DMF",  "OT_SPF",   "OT_SPF",  "OT_SPF",   "OT_SPF",   "SDN_DEM",
     "SDN_DEM",   "SDN_DEM", "SDN_DEM",  "SSC_DEM",   "SSC_DEM",  "TBB_CRU",   "TBB_DMF")  
   tacsatp$LE_MET         <- as.character(tacsatp$LE_MET)
   } else{
    if(a_year=="2011"){
     levels(tacsatp$LE_MET) <-   c(  ## REPLACE LEVELS WITH CAUTION ## adapt to your own list!!
      "DRB_MOL",      "NA",   "OT_CRU", "OT_CRU",   "OT_CRU",  "OT_SPF", "OT_DMF", "OT_SPF",  
      "OT_SPF",   "OT_DMF",  "OT_MIX_NEP",   "OT_MIX_NEP", "OT_MIX_NEP",   "OT_MIX_NEP",  "OT_SPF",  "OT_SPF",  
      "OT_SPF",   "OT_SPF",   "OT_SPF",   "OT_CRU",  "OT_SPF",     "OT_DMF", "OT_SPF",   "OT_MIX_NEP",  
      "OT_MIX_NEP",  "OT_SPF",   "OT_SPF",  "OT_SPF",   "OT_SPF",   "SDN_DEM", "SDN_DEM",   "SDN_DEM",
      "SDN_DEM",  "SSC_DEM", "SSC_DEM",   "TBB_CRU",   "TBB_DMF",     "TBB_DMF")  
    }else{
     if(a_year=="2012") {
       levels(tacsatp$LE_MET) <-   c("DRB_MOL", "NA", "OT_CRU", "OT_CRU", "OT_CRU",  "OT_DMF", "OT_DMF", "OT_DMF",
        "OT_DMF", "OT_DMF", "OT_MIX_NEP", "OT_MIX_NEP",  "OT_MIX_NEP", "OT_MIX_NEP",
        "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF", "OT_DMF", "OT_MIX_NEP",
        "OT_MIX_NEP",  "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF", 
        "SDN_DEM", "SDN_DEM",   "SDN_DEM", "SDN_DEM", "SSC_DEM", "SSC_DEM", "TBB_CRU", "TBB_DMF")
    } else{
      if(a_year=="2013") {
       levels(tacsatp$LE_MET) <-   c("DRB_MOL", "NA", "OT_CRU", "OT_DMF", "OT_DMF", "OT_DMF", "OT_DMF", "OT_DMF", 
       "OT_MIX_NEP", "OT_MIX_NEP", "OT_MIX_NEP", "OT_MIX_NEP", "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF", "OT_DMF",
       "OT_DMF", "OT_MIX_NEP", "OT_SPF", "OT_SPF", "OT_DMF",  "OT_SPF",  "SDN_DEM", "SDN_DEM", "SDN_DEM", "SDN_DEM",
       "SSC_DEM", "SSC_DEM",   "TBB_CRU")
       
    } else{
    stop('adapt the BENTHIS metiers for this year')
    }
   }
  }  
  }
  }
  initVersusBenthisMetiers <-  tacsatp [!duplicated(data.frame(tacsatp$LE_MET_init, tacsatp$LE_MET)), 
                                    c('LE_MET_init', 'LE_MET')]
  save(initVersusBenthisMetiers, file=file.path(outPath,a_year,"initVersusBenthisMetiers.RData"))



  
  #-----------------------------------------------------------------------------
  # Define activity
  #-----------------------------------------------------------------------------
  
  idx               <- which(is.na(tacsatp$VE_REF) == T   | is.na(tacsatp$SI_LONG) == T | is.na(tacsatp$SI_LATI) == T |
                             is.na(tacsatp$SI_DATIM) == T |  is.na(tacsatp$SI_SP) == T)
  if(length(idx)>0) tacsatp         <- tacsatp[-idx,]

  if(.Platform$OS.type == "windows" && TRUE) {
    storeScheme       <- activityTacsatAnalyse(tacsatp, units = "year", analyse.by = "LE_GEAR",identify="means")
    storeScheme       <- storeScheme[which(is.na(storeScheme$analyse.by)==F),]

    storeScheme$years <- as.numeric(as.character(storeScheme$years))
    storeScheme       <- storeScheme[storeScheme$years==a_year,]
    save(storeScheme, file=file.path(outPath,a_year,"storeScheme.RData"))
  }  else{
    load(file.path(outPath,a_year,"storeScheme.RData"))
  }

  tacsatp$year      <- format(tacsatp$SI_DATIM, "%Y")
  require(mixtools)
  activity          <- activityTacsat(tacsatp,units="year",analyse.by="LE_GEAR",storeScheme,
                              plot=FALSE, level="all")
  tacsatp$SI_STATE  <- NA
  tacsatp$SI_STATE  <- activity

  #- Plot the result
  if(FALSE){
    result          <- table(tacsatp$SI_STATE,tacsatp$SI_SP,tacsatp$LE_GEAR)
    par(mfrow=rep(ceiling(sqrt(length(unique(tacsatp$LE_GEAR)))),2))
    for(i in 1:dim(result)[3])
      barplot(result[,,i],col=1:3)
  }

  #- General speed rules for remaining records
  idx               <- which(is.na(tacsatp$SI_STATE))
  if(length(idx)>0){
    tacsatp$SI_STATE[idx[which(tacsatp$SI_SP[idx] >= 1.5 &
                               tacsatp$SI_SP[idx] <= 7.5)]] <- 'f'
    tacsatp$SI_STATE[idx[which(tacsatp$SI_SP[idx] <  1.5)]] <- 'h'
    tacsatp$SI_STATE[idx[which(tacsatp$SI_SP[idx] >  7.5)]] <- 's'
  }
  save(tacsatp,     file=file.path(outPath,a_year,"tacsatActivity.RData"))

  # Labelling each haul (caution: to do before discarding the steaming points...)
  tacsatp   <- labellingHauls(tacsatp)


  #-----------------------------------------------------------------------------
  # Interpolation (of fishing sequences only)
  #-----------------------------------------------------------------------------
  dir.create(file.path(outPath,a_year,"interpolated"))
  tacsatp           <- orderBy(~VE_REF+SI_DATIM,data=tacsatp)

  # KEEP ONLY fish. seq. bounded by steaming points
  tacsatp$SI_STATE_num <- NA
  tacsatp$SI_STATE_num[which(tacsatp$SI_STATE=="h")] <- 1
  tacsatp$SI_STATE_num[tacsatp$SI_STATE=="f"] <- 2
  tacsatp$SI_STATE_num[tacsatp$SI_STATE=="s"] <- 3
  is_transition     <- c(0,diff(tacsatp$SI_STATE_num))
  is_transition2    <- c(diff(tacsatp$SI_STATE_num), 0)
  tacsatp           <- tacsatp[ !is.na(tacsatp$SI_STATE_num) & (tacsatp$SI_STATE_num ==2 |
                                      is_transition!=0 | is_transition2!=0),]
  tacsatp           <- tacsatp[,-grep("SI_STATE_num",colnames(tacsatp))]
  tacsatp$SI_STATE  <- "f"
  
  #- Gear specific fm parameters (parameters tuned with st=c(4,8))
  #- (should be informed for all below towedGears)
  fm        <- list(TBB=0.4575222,
                    OTB=0.4217848,
                    OTT=0.4217848,
                    PTB=0.4217848,
                    HMD=0,
                    DRB=0)

  # per gear per vessel
  for(iGr in towedGears){
    tacsatpGear             <- tacsatp[!is.na(tacsatp$LE_GEAR) & tacsatp$LE_GEAR==iGr,]

    for(iVE_REF in sort(unique(tacsatpGear$VE_REF))){
      tacsatpGearVEREF      <- tacsatpGear[tacsatpGear$VE_REF %in% iVE_REF,]
      if(nrow(tacsatpGearVEREF)>3){
        cat(paste(iGr, " ", iVE_REF, "\n"))

        #Interpolate according to the cubic-hermite spline interpolation
        interpolationcHs    <- interpolateTacsat(tacsatpGearVEREF,
                                                 interval= VMS_ping_rate_in_hour*60,   # THE PING RATE IS COUNTRY-SPECIFIC ##
                                                 margin  = round(VMS_ping_rate_in_hour*60*0.1), # i.e. will make disconnected interpolations if interval out of the 50 70min range
                                                 res     = round(VMS_ping_rate_in_hour*(100/2)),
                                                 method  = "cHs",
                                                 params  = list(fm=fm[[iGr]],distscale=0,sigline=0,st=c(4,8)),   # rmenber that st not in use....
                                                 headingAdjustment=0,
                                                 fast    = FALSE)

        # Get the ranges of the total picture
        if(FALSE){
          ranges <- do.call(rbind,lapply(interpolationcHs,function(x){return(apply(x[-1,],2,range))}))
          xrange <- range(ranges[,"x"])
          yrange <- range(ranges[,"y"])

          plot(tacsatpGearVEREF$SI_LONG, tacsatpGearVEREF$SI_LATI,
               xlim=xrange,ylim=yrange,pch=19,col="blue",xlab="Longitude",ylab="Latitude")
          for(iInt in 1:length(interpolationcHs))
            lines(interpolationcHs[[iInt]][-1,1],interpolationcHs[[iInt]][-1,2])
        }

        # Convert the interpolation to tacsat style data
        medx               <- median(tacsatpGearVEREF$SI_LONG,na.rm=T); medy <- median(tacsatpGearVEREF$SI_LATI,na.rm=T)
        npoints            <- ceiling(median(tacsatpGearVEREF$SI_SP,na.rm=T) * 1.852 * VMS_ping_rate_in_hour /
                                      mean(distance(medx,medy,medx+1/60,medy),distance(medx,medy,medx,medy+1/60))) + 1
        tacsatIntGearVEREF <- interpolation2Tacsat(interpolationcHs, tacsatpGearVEREF,npoints=ifelse(npoints<2,2,npoints))
        tacsatIntGearVEREF <- tacsatIntGearVEREF[!duplicated(apply(tacsatIntGearVEREF[,c("SI_LONG","SI_LATI","SI_DATIM")],1,paste,collapse="_")),]

  
        save(tacsatIntGearVEREF, file=file.path(outPath,a_year,"interpolated",
                                                paste("tacsatSweptArea_",iVE_REF, "_", iGr, ".RData", sep="")),compress=T)
      }
    }
  }


  # per gear per vessel
  for(iGr in seineGears){
    tacsatpGear        <- tacsatp[!is.na(tacsatp$LE_GEAR) & tacsatp$LE_GEAR==iGr,]

    for(iVE_REF in sort(unique(tacsatpGear$VE_REF))){
      cat(paste(iGr, " ", iVE_REF, "\n"))
      tacsatpGearVEREF <- tacsatpGear[tacsatpGear$VE_REF %in% iVE_REF,]
      tacsatpGearVEREF <- tacsatpGearVEREF[tacsatpGearVEREF$SI_STATE=='f',] # keep fishing pings only

  
      tacsatIntGearVEREF <- tacsatpGearVEREF

      save(tacsatIntGearVEREF, file=file.path(outPath,a_year,"interpolated",
                                              paste("tacsatSweptArea_",iVE_REF, "_", iGr, ".RData", sep="")))
    }
  }
} # end TRUE/FALSE

 
 
 #-----------------------------------------------------------------------------
 # compute (discrete point) effort_days and effort_KWdays
 #-----------------------------------------------------------------------------

  if(FALSE){
  library(doBy)
  tacsatp                  <- orderBy(~VE_REF+SI_DATIM+FT_REF,data=tacsatp)
  tacsatp$effort_days      <- as.numeric(as.character(difftime(c(tacsatp$SI_DATIM[-1],0),tacsatp$SI_DATIM,units="days")))
  tacsatp$effort_KWdays    <- tacsatp$effort_days *  as.numeric(as.character(tacsatp$VE_KW))
  #tacsatp$effort_days   [which(diff(as.numeric(tacsatp$HL_ID))<0)]     <- 0  # correct (to do: but pble of steaming pts which have been converted to fish. pts)
  #tacsatp$effort_KWdays [which(diff(as.numeric(tacsatp$HL_ID))<0)]     <- 0  # correct
  ## for the time being, do:
  tacsatp$effort_days[tacsatp$effort_days>0.014] <- 0  # correct (i.e. set at 0 if >3hours as a sign for a change of haul)
  }

  #-----------------------------------------------------------------------------
  # Add "gear width-vessel size" relationships table of parameters.
  #-----------------------------------------------------------------------------
  #gear_param_per_metier       <- read.table(file=file.path(dataPath, "estimates_for_gear_param_per_metier.txt"))
  # an equivalent is:
  
  gear_param_per_metier <- data.frame(
  a_metier=c('OT_CRU','OT_CRU','OT_DMF','OT_DMF','OT_MIX','OT_MIX','OT_MIX_ARA','OT_MIX_ARA','OT_MIX_DMF_BEN','OT_MIX_DMF_BEN','OT_MIX_DMF_PEL','OT_MIX_DMF_PEL','OT_MIX_DPS','OT_MIX_DPS','OT_MIX_NEP','OT_MIX_NEP','OT_MIX_TGS_CTC','OT_MIX_TGS_CTC','OT_MIX_TGS_OCC','OT_MIX_TGS_OCC','OT_SPF','OT_SPF','TBB_CRU','TBB_CRU','TBB_DMF','TBB_DMF','TBB_MOL','TBB_MOL','DRB_MOL','DRB_MOL','SDN_DEM','SDN_DEM','SSC_DEM','SSC_DEM'),
  param=c('a','b','a','b','a','b','a','b','a','b','a','b','a','b','a','b','a','b','a','b','a','b','a','b','a','b','a','b','a','b','a','b','a','b'),
  Estimate=c(5.10393560454806,0.468985756915913,9.6053549509854,0.433672763959314,10.6607888271164,0.292055014993337,37.5271604597435,0.149004797319136,3.21410379943408,77.981158829069,6.63707197355847,0.770594580782091,26.6738247840508,0.210221545999405,3.92727763464472,35.8253721834011,6.23686411376723,0.767375050454527,0.0192465419797634,119.140335982507,0.965238378524667,68.3889717127507,1.48117115311386,0.457788539321641,0.660086393453441,0.507845311175148,0.953001905566232,0.709356826689359,0.314245137194503,1.24544036138755,1948.83466676682,0.236271746198865,4461.27004311913,0.117589220782479),
  Std..Error=c(1.81527145191998,0.0597519960969362,3.98228885098937,0.067572002767068,6.69386377505425,0.104413257104915,10.6717875588847,0.044963446750424,1.67854244656697,40.9297885227685,2.69086696344053,0.126123213329976,5.37466576335144,0.030829495804396,0.928442484509969,21.0228522096513,1.46159830273852,0.0732116002636393,0.000552819642352548,0.510207569180525,0.205245990518183,7.45180177818494,0.278399892100703,0.0346555048025894,0.172902115850281,0.0388684340513048,0.315715856194751,0.138412196798781,0.110027479611801,0.10614681568516,637.25152416296,0.0636712369543136,1665.50234108383,0.118756519107319),
  t.value=c(2.81166521907769,7.84887179593252,2.41201864314765,6.41793562718951,1.59262112068153,2.79710664230959,3.51648308708138,3.31390958851994,1.91481830322951,1.90524216331315,2.46651806415295,6.10985527910701,4.96288066244663,6.81884476260001,4.22996329893018,1.70411568450042,4.26715336360309,10.4816046595234,34.8152281598731,233.513462322532,4.70283670871103,9.17750817164227,5.32030074414718,13.2096918492278,3.81768835047121,13.0657517744299,3.01854305657162,5.12495894939517,2.8560604887363,11.733186279291,3.05818753329251,3.71080816866175,2.67863330664306,0.990170658978435),
  Pr...t..=c(0.00613312535554725,1.21619365805854e-11,0.021410083292817,2.48114253493853e-07,0.114790848188445,0.00631861326022122,0.000513087659147687,0.0010462790834138,0.0692370736030276,0.0705334706657513,0.0147045751318625,7.39218704723967e-09,1.2637878625965e-05,2.97113026239585e-08,0.000166717383514359,0.097483711710908,0.000314181622785133,5.0948672020349e-10,9.05842416252619e-12,5.10054218622276e-20,0.000204968683311441,5.36482029322678e-08,0.00313939649832079,4.44157761915604e-05,0.000458495488420268,5.11509704563588e-16,0.00678642704689924,5.16047183433098e-05,0.0075895814688592,6.18091407283774e-13,0.00391206507124518,0.000614325243514857,0.0438919330122769,0.367557330382699),
  equ=c('DoS=a*(kW^b)','DoS=a*(kW^b)','DoS=a*(kW^b)','DoS=a*(kW^b)','DoS=a*(kW^b)','DoS=a*(kW^b)','DoS=a*(kW^b)','DoS=a*(kW^b)','DoS=(a*LOA)+b','DoS=(a*LOA)+b','DoS=a*(LOA^b)','DoS=a*(LOA^b)','DoS=a*(kW^b)','DoS=a*(kW^b)','DoS=(a*LOA)+b','DoS=(a*LOA)+b','DoS=a*(LOA^b)','DoS=a*(LOA^b)','DoS=(a*kW)+b','DoS=(a*kW)+b','DoS=(a*LOA)+b','DoS=(a*LOA)+b','beamw=a*(kW^b)','beamw=a*(kW^b)','beamw=a*(kW^b)','beamw=a*(kW^b)','beamw=a*(LOA^b)','beamw=a*(LOA^b)','dredgew=a*(LOA^b)','dredgew=a*(LOA^b)','seineropel=a*(kW^b)','seineropel=a*(kW^b)','seineropel=a*(LOA^b)','seineropel=a*(LOA^b)'),
  nb_records=c(124,124,39,39,94,94,271,271,48,48,190,190,45,45,53,53,24,24,12,12,19,19,7,7,42,42,22,22,33,33,47,47,8,8)
  )
  
 
 
  


#-----------------------------------------------------------------------------
# Create one swept area dataset
#-----------------------------------------------------------------------------
fls <- dir(file.path(outPath,a_year,"interpolated"))

lst <- list(); count <- 0
cols2keep <- c("SI_LATI","SI_LONG","SI_DATE","LE_GEAR","LE_MET","SWEPT_AREA_KM2","SWEPT_AREA_KM2_LOWER","SWEPT_AREA_KM2_UPPER")
for(iFile in fls){
  cat(paste(iFile, "\n"))
  count <- count+1
  load(file.path(outPath,a_year,"interpolated",iFile))


  
  #- Make selection for gears where you already have gear width and which not
  ctry <- "XXX"
  if(ctry=="NLD"){
  # compute the swept area
  tacsatIntGearVEREF <- compute_swept_area (tacsatIntGearVEREF, gear_param_per_metier, towedGears, seineGears, VMS_ping_rate_in_hour, already_informed_width_for=c('DRB', 'TBB'))
  } else{
  # compute the swept area
  tacsatIntGearVEREF <- compute_swept_area (tacsatIntGearVEREF, gear_param_per_metier, towedGears, seineGears, VMS_ping_rate_in_hour, already_informed_width_for=NULL)
  }
  
  if(nrow(tacsatIntGearVEREF[is.na(tacsatIntGearVEREF$SWEPT_AREA_KM2),])!=0 ) print('check for NAs')
  # NAs are acceptable if the metier was not informed (e.g. No_Matrix6)...
  print(unique(tacsatIntGearVEREF[is.na(tacsatIntGearVEREF$SWEPT_AREA_KM2),"LE_MET_init"]))
 
  lst[[count]] <- tacsatIntGearVEREF[,cols2keep]

}
tacsatSweptArea   <- do.call(rbind,lst)

# check NAs (approx. 2% of the records)
nrow(tacsatSweptArea[is.na(tacsatSweptArea$SWEPT_AREA_KM2),])

# save
save(tacsatSweptArea, file=file.path(outPath,a_year, paste("tacsatSweptArea.RData", sep="")),compress=T)





#-----------------------------------------------------------------------------
# GRIDDING (IN DECIMAL DEGREES OR UTM COORD)
#-----------------------------------------------------------------------------
  # using a quick gridding code at various resolution.

  
  # For example, grid the swept area, or the number of hauls, or the fishing pressure, etc.
  sh1 <<- readShapePoly(file.path(polPath,"francois_EU"))


  gridding <- function(tacsatp, 
                      general=list(what="SWEPT_AREA_KM2", 
                                   a_func="sum",      #a_func  <- function(x) {unique(length(x))} # for HL_ID, to be tested.
                                   is_utm=TRUE,
                                   utm_zone=32,
                                   BENTHISmetiers=c('DRB_MOL', 'NA', 'OT_CRU', 'OT_SPF', 'OT_DMF', 'OT_MIX_NEP', 'SDN_DEM', 'SSC_DEM', 'TBB_CRU', 'TBB_DMF'),
                                   we=10,
                                   ea=15,
                                   so=52.5,
                                   no=59,
                                   the_breaks =c(0, (1:12)^2.5 ),
                                   output_file_name='GriddedSweepAreaExample',
                                   export_cumul=FALSE # because time consuming....
                               )
                       ){
                   
    # subset for relevant fisheries
    this            <- tacsatp [tacsatp$LE_MET %in% general$BENTHISmetiers, ]

    # restrict the study area
    this <- this[this$SI_LONG>general$we & this$SI_LONG<general$ea & this$SI_LATI>general$so & this$SI_LATI<general$no,]

    # grid the data (in decimal or in UTM)
    if(general$is_utm){
      dx <- 0.0002 # 5km 
      dx <- 0.001 # 1km 
      # convert to UTM
      library(sp)
      library(rgdal)
      SP <- SpatialPoints(cbind(as.numeric(as.character(this$SI_LONG)), as.numeric(as.character(this$SI_LATI))),
                       proj4string=CRS("+proj=longlat +datum=WGS84"))
      this <- cbind(this,
                 spTransform(SP, CRS(paste("+proj=utm  +ellps=intl +zone=",utm_zone," +towgs84=-84,-107,-120,0,0,0,0,0", sep=''))))    # convert to UTM
      this            <- this [, c('SI_LONG', 'SI_LATI', 'SI_DATE', 'coords.x1', 'coords.x2', what)]
      this$round_long <- round(as.numeric(as.character(this$coords.x1))*dx)
      this$round_lat  <- round(as.numeric(as.character(this$coords.x2))*dx)
      this            <- this[, !colnames(this) %in% c('coords.x1', 'coords.x2')]
      this$cell       <- paste("C_",this$round_long,"_", this$round_lat, sep='')
      this$xs         <- (this$round_long/(dx))
      this$ys         <- (this$round_lat/(dx))
      colnames(this) <- c('x', 'y', 'date', 'what', 'round_long', 'round_lat', 'cell', 'xs', 'ys')  
      # retrieve the geo resolution, for info
      # unique(this$ys)
      # unique(this$xs)
    }  else {
      dx <- 20
      this <- this [, c('SI_LONG', 'SI_LATI', 'SI_DATE', what)]
      this$round_long <- round(as.numeric(as.character(this$SI_LONG))*dx*2)
      this$round_lat  <- round(as.numeric(as.character(this$SI_LATI))*dx)
      this$cell       <- paste("C_",this$round_long,"_", this$round_lat, sep='')
      this$xs         <- (this$round_long/(dx*2))
      this$ys         <- (this$round_lat/(dx))
      colnames(this) <- c('x', 'y', 'date', 'what', 'round_long', 'round_lat', 'cell', 'xs', 'ys')  
      # retrieve the geo resolution in degree, for info
      #long <- seq(1,15,by=0.01)
      #res_long <- diff( long [1+which(diff(round(long*dx)/dx/2)!=0)] )
      #res_lat <- diff( long [1+which(diff(round(long*dx/2)/dx)!=0)] )
      #print(res_long) ; print(res_lat)
      }
     # if the coordinates in decimal then dx=20 corresponds to grid resolution of 0.05 degrees
     # i.e. a 3´ angle = 3nm in latitude but vary in longitude (note that a finer grid will be produced if a higher value for dx is put here)
     # if coord in UTM then 0.001 correspond to grid of 1 by 1 km (to check)



    # a quick gridding method...
    background <- expand.grid(
                              x=0,
                              y=0,
                              date=0,  
                              what=0,
                              round_long=seq(range(this$round_long)[1], range(this$round_long)[2], by=1),
                              round_lat=seq(range(this$round_lat)[1], range(this$round_lat)[2], by=1),
                              cell=0,
                              xs=0,
                              ys=0
                              )
    this <- rbind(this, background)
    the_points <- tapply(this$what,
                  list(this$round_lat, this$round_long), general$a_func, na.rm=TRUE)

    if(!is_utm) xs <- (as.numeric(as.character(colnames(the_points)))/(dx*2))
    if(is_utm) xs <- (as.numeric(as.character(colnames(the_points)))/(dx))
    ys <- (as.numeric(as.character(rownames(the_points)))/(dx))

   
    graphics:::image(
     x=xs,
     y=ys,
     z= t(the_points)  ,
     breaks=c(general$the_breaks),
     col = terrain.colors(length(general$the_breaks)-1),
     useRaster=FALSE,
     xlab="",
     ylab="",
     axes=FALSE,
     xlim=range(xs), ylim=range(ys),
     add=FALSE
     )
    title("")

    # land
    sh1 <- readShapePoly(file.path(polPath,"francois_EU"),  proj4string=CRS("+proj=longlat +datum=WGS84"))
    if(is_utm) sh1 <- spTransform(sh1, CRS(paste("+proj=utm  +ellps=intl +zone=",general$utm_zone," +towgs84=-84,-107,-120,0,0,0,0,0", sep='')))
    plot(sh1, add=TRUE, col=grey(0.7))

    legend("topright", fill=terrain.colors(length(general$the_breaks)-1), title=paste(what, " 1 x 1km cell"),
             legend=round(general$the_breaks[-1],1), bty="n", cex=0.8, ncol=2)
    box()
    axis(1)
    axis(2, las=2)

    if(general$is_utm){
      mtext(side=1, "Eastings", cex=1, adj=0.5, line=2)
      mtext(side=2, "Northings", cex=1, adj=0.5, line=2)
    } else{
      mtext(side=1, "Longitude", cex=1, adj=0.5, line=2)
      mtext(side=2, "Latitude", cex=1, adj=0.5, line=2)
      #points (tacsatp [,c('SI_LONG', 'SI_LATI')], pch=".", col="white")
    }


    # save
    savePlot(filename=file.path(outPath,a_year, paste(general$output_file_name,".jpeg", sep='')), type="jpeg")
 
  

    # export the quantity per cell and date
    library(data.table)
    DT                <-  data.table(this)
    qu                <-  quote(list(sum(what)))
    quantity_per_cell <- DT[,eval(qu), by=list(cell, xs,ys)]
    quantity_per_date <- DT[,eval(qu), by=list(date, xs,ys)]
    quantity_per_cell_date <- DT[,eval(qu), by=list(cell,date, xs,ys)]
    quantity_per_cell_date <- as.data.frame(quantity_per_cell_date)
    colnames(quantity_per_cell_date)   <- c('cell','date', 'xs','ys', 'quantity')
    rm(DT) ; gc(reset=TRUE)
    save(quantity_per_cell_date, res_long, res_lat,  general, file=file.path(outPath,a_year, paste(general$output_file_name,"_quantity_per_cell_date.RData", sep='')) )


    # export the cumul per cell and date
    if(general$export_cumul){
    quantity_per_cell_date <- orderBy(~date, data=quantity_per_cell_date)
    quantity_cumul_per_cell_date <- do.call("rbind", lapply(
      split(quantity_per_cell_date, f=quantity_per_cell_date$cell),
               function(x){
               x$quantity <- cumsum(x$quantity)
               x
               })  )
    # check the cumul on a given cell
    # head(quantity_cumul_per_cell_date[quantity_cumul_per_cell_date$cell=="C_197_2722",])               
    save(quantity_cumul_per_cell_date, res_long, res_lat,  general, file=file.path(outPath,a_year, paste(general$output_file_name,"_quantity_cumul_per_cell_date.RData",sep='')) )
    }

  cat(paste("Please find the saved map and export files in ",file.path(outPath,a_year), "\n"))
 return() 
 }
 
 #-------------------------------------
 # CALLS for mapping and exporting-----
  load(file.path(outPath, a_year, "tacsatSweptArea.RData")) # get 'tacsatSweptArea' with all data
  tacsatp <- tacsatSweptArea 
  graphics.off()
  
  gridding(tacsatp=tacsatp,   general=list(what="SWEPT_AREA_KM2", 
                                   a_func="sum",      #a_func  <- function(x) {unique(length(x))} # for HL_ID, to be tested.
                                   is_utm=TRUE,
                                   utm_zone=32,
                                   BENTHISmetiers=c('DRB_MOL', 'NA', 'OT_CRU', 'OT_SPF', 'OT_DMF', 'OT_MIX_NEP', 'SDN_DEM', 'SSC_DEM', 'TBB_CRU', 'TBB_DMF'),
                                   we=10,
                                   ea=15,
                                   so=52.5,
                                   no=59,
                                   the_breaks =c(0, (1:12)^1.5 ),
                                   output_file_name='GriddedSweepAreaExample_allmet',
                                   export_cumul=FALSE
                               ))
 
 
   gridding(tacsatp=tacsatp,   general=list(what="SWEPT_AREA_KM2", 
                                   a_func="sum",      #a_func  <- function(x) {unique(length(x))} # for HL_ID, to be tested.
                                   is_utm=TRUE,
                                   utm_zone=32,
                                   BENTHISmetiers=c('SDN_DEM', 'SSC_DEM'),
                                   we=10,
                                   ea=15,
                                   so=52.5,
                                   no=59,
                                   the_breaks =c(0, (1:12)^2.5 ),
                                   output_file_name='GriddedSweepAreaExample_seinersonly',
                                   export_cumul=FALSE
                               ))
 
 
   gridding(tacsatp=tacsatp,   general=list(what="SWEPT_AREA_KM2", 
                                   a_func="sum",      #a_func  <- function(x) {unique(length(x))} # for HL_ID, to be tested.
                                   is_utm=TRUE,
                                   utm_zone=32,
                                   BENTHISmetiers=c('DRB_MOL', 'NA', 'OT_CRU', 'OT_SPF', 'OT_DMF', 'OT_MIX_NEP', 'SDN_DEM', 'SSC_DEM', 'TBB_CRU', 'TBB_DMF'),
                                   we=1,
                                   ea=15,
                                   so=52.5,
                                   no=63,
                                   the_breaks =c(0, (1:12)^1.5 ),
                                   output_file_name='GriddedSweepAreaExample_allmet_largerarea',
                                   export_cumul=FALSE
                               ))



 #-----------------------------------------------------------------------------
 # LINK TO HABITAT MAPS
 #-----------------------------------------------------------------------------


  #------------------
  # SHAPEFILE--------
  #------------------
  library(maptools)
  library(rgdal)
  
  # load a habitat map shape file
  habitat_map           <- readShapePoly(file.path(polPath,"sediment_lat_long"),
                                         proj4string=CRS("+proj=longlat +ellps=WGS84"))

  sh_coastlines            <- readShapePoly(file.path(polPath,"francois_EU"))


  load(file.path(outPath, a_year, "tacsatSweptArea.RData")) # get 'tacsatSweptArea' with all data
  #....or load only one instance eg load(file.path(outPath, a_year, "interpolated","tacsatSweptArea_DNK000001744_OTB.RData"))
  #tacsatSweptArea <- tacsatIntGearVEREF

  coord <-  tacsatSweptArea[, c('SI_LONG', 'SI_LATI')]
  names(habitat_map) # return the name of the coding variable

  #Turn the polygons into spatial polygons
  sp <- SpatialPolygons(habitat_map@polygons, proj4string=CRS("+proj=longlat +ellps=WGS84"))
 
  #Turn the point into a spatial point
  spo <- SpatialPoints(data.frame(SI_LONG=coord[,1], SI_LATI=coord[,2]), proj4string=CRS("+proj=longlat +ellps=WGS84"))         
 
  #Use the magic 'over' function to see in which polygon it is located
  idx <- over(spo,sp); print(idx)
  tacsatSweptArea$SUBSTRATE <- habitat_map$BAL_CODE[idx]

  # plot
  plot(habitat_map, xlim=c(11,14), ylim=c(55,56))
  axis(1) ; axis(2, las=2) ; box()
  points(tacsatSweptArea[, c("SI_LONG","SI_LATI")], col=tacsatSweptArea$SUBSTRATE, pch=".")
  #plot(sh_coastlines, add=TRUE)

  # save
  savePlot(filename=file.path(outPath,a_year,"VMSpingsAttachedToSedimentMap.jpeg"), type="jpeg")


  #------------------
  # RASTER-----------
  #------------------
  sh_coastlines            <- readShapePoly(file.path(polPath,"francois_EU"))

  ## use point-raster overlay.......
  library(raster)
  landscapes       <- raster(file.path(polPath, "landscapes.tif"))    # probably need an update of rgdal here....
  newproj          <- "+proj=longlat +datum=WGS84"
  landscapes_proj  <- projectRaster(landscapes, crs=newproj)

  save(landscapes_proj, file=file.path(polPath,a_year, "landscapes_proj.RData"))

  load(file.path(polPath, "landscapes_proj.RData"))

  load(file.path(outPath, a_year, "tacsatSweptArea.RData")) # get 'tacsatp' with all data
  #....or load only one instance eg load("C:\\merging\\BENTHIS\\outputs\\interpolated\\tacsatSweptArea_DNK000005269_OTB.RData"))
  #tacsatp <- tacsatpGearVEREF


  coord <- cbind(x=anf(tacsatp$SI_LONG), y=anf(tacsatp$SI_LATI))

  dd <- extract (landscapes_proj, coord[,1:2]) # get the landscape on the coord points!

  coord <- cbind(coord,  landscapes_code=cut(dd, breaks=c(0,100,200,300,400,500,600)))

  tacsatp <- cbind(tacsatp, landscapes_code= coord[,'landscapes_code'])

  # plot and save...
  plot(landscapes_proj, xlim=c(10,14), ylim=c(54.5,56.5))
  plot(sh_coastlines,  xlim=c(10,14), ylim=c(54.5,56.5), add=TRUE)  # better for plotting the western baltic sea coastline!
  points(coord[,"x"], coord[,"y"], col=coord[,"landscapes_code"], pch=".", cex=1)

  # save
  savePlot(filename=file.path(outPath,a_year, "VMSpingsAttachedToLandscapeMap.jpeg"), type="jpeg")



 #-----------------------------------------------------------------------------
 #  SEVERITY OF IMPACT
 #-----------------------------------------------------------------------------

  # create a fake input file to show the required format
  dd <- tacsatp[!duplicated(data.frame(tacsatp$VE_REF,tacsatp$LE_GEAR, tacsatp$LE_MET)), ]
  fake_gear_metier_habitat_severity_table <- dd[,c('VE_REF', 'LE_GEAR', 'LE_MET')]
  fake_gear_metier_habitat_severity_table <- cbind(fake_gear_metier_habitat_severity_table, HAB_SEVERITY=1)
  gear_metier_habitat_severity_table  <- fake_gear_metier_habitat_severity_table
  gear_metier_habitat_severity_table  <-   gear_metier_habitat_severity_table[complete.cases( gear_metier_habitat_severity_table),]
  save(gear_metier_habitat_severity_table,   file=paste(dataPath,a_year,"gear_metier_habitat_severity_table.RData",   sep=""))

  # load a table for HAB_SEVERITY per vid LE_REF per gr LE_GEAR per met LE_MET
  load(file.path(dataPath, "gear_metier_habitat_severity_table.RData"))
  tacsatp <- merge(tacsatp, gear_metier_habitat_severity_table)
  save(tacsatp,   file=paste(outPath,a_year,"tacsatMergedHabSeverity.RData",   sep=""))

  # pressure: weigh the swept area by the severity
  tacsatp$pressure <- tacsatp$HAB_SEVERITY * tacsatp$SWEPT_AREA_KM2


























