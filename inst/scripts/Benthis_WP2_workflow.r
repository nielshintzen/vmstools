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


#a_year      <- 2005
#a_year      <- 2006
#a_year      <- 2007
#a_year      <- 2008
#a_year      <- 2009
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
        # correct the transition between sequential fishing events
        idx <- which(diff(tacsatIntGearVEREF$SI_DATIM)/60 > 15)   # if interval > 15 min then points belong to a different fishing event
        tacsatIntGearVEREF[ idx, c('SWEPT_AREA_KM2', 'SWEPT_AREA_KM2_LOWER', 'SWEPT_AREA_KM2_UPPER')] <- NA
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
     if(a_year=="2005"){
     levels(tacsatp$LE_MET) <-   c(  ## REPLACE LEVELS WITH CAUTION ## adapt to your own list!!
     "DRB_MOL", "OT_DMF", "NA", "OT_CRU",  "OT_CRU",  "OT_CRU", "OT_CRU",  "OT_CRU", "OT_CRU",  "OT_CRU",
     "OT_MIX_DMF_PEL", "OT_DMF", "OT_DMF", "OT_DMF", "OT_MIX_DMF_PEL",   "OT_MIX_DMF_PEL",  "OT_DMF", "OT_DMF",  "OT_DMF", "OT_DMF", 
     "OT_SPF",  "OT_SPF",  "OT_SPF",  "OT_SPF",  "OT_SPF",  "OT_CRU",  "OT_MIX_DMF_PEL",  "OT_DMF", "OT_DMF", "OT_DMF",
     "OT_MIX_DMF_PEL",  "OT_DMF", "OT_SPF",  "OT_SPF",  "OT_SPF", "OT_SPF",  "OT_SPF",   "SDN_DEM", "SDN_DEM",  "SDN_DEM",
     "SDN_DEM",  "SSC_DEM",   "SSC_DEM","SSC_DEM",   "TBB_CRU",   "TBB_DMF", "TBB_DMF")  
   } else{ 
      if(a_year=="2006"){
     levels(tacsatp$LE_MET) <-   c(  ## REPLACE LEVELS WITH CAUTION ## adapt to your own list!!
      "DRB_MOL", "NA",  "OT_CRU", "OT_CRU",    "OT_CRU", "OT_CRU", "OT_CRU","OT_CRU","OT_CRU",  "OT_CRU", 
      "OT_MIX_DMF_PEL",     "OT_DMF", "OT_DMF",   "OT_DMF", "OT_MIX_DMF_PEL", "OT_MIX_DMF_PEL",  "OT_DMF","OT_DMF", "OT_DMF", "OT_DMF", 
      "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF",  "OT_CRU",  "OT_CRU", "OT_MIX_DMF_PEL", "OT_DMF",
      "OT_DMF", "OT_DMF", "OT_MIX_DMF_PEL", "OT_MIX_DMF_PEL",  "OT_DMF",  "OT_SPF","OT_SPF","OT_SPF", "OT_SPF", "OT_SPF",  
      "OT_SPF",  "SDN_DEM", "SDN_DEM", "SDN_DEM", "SDN_DEM",  "SSC_DEM", "SSC_DEM", "SSC_DEM", "SSC_DEM", "SSC_DEM", 
      "TBB_CRU",   "TBB_DMF",   "TBB_DMF", "TBB_DMF")  

   } else{ 
   if(a_year=="2007"){
     levels(tacsatp$LE_MET) <-   c(  ## REPLACE LEVELS WITH CAUTION ## adapt to your own list!!
     "DRB_MOL", "NA", "OT_CRU", "OT_CRU", "OT_CRU", "OT_CRU","OT_CRU","OT_CRU",  "OT_MIX_DMF_PEL",  "OT_DMF",
     "OT_DMF",   "OT_DMF", "OT_MIX_DMF_PEL",   "OT_MIX_DMF_PEL",   "OT_DMF",   "OT_DMF",  "OT_DMF",  "OT_SPF", "OT_SPF", "OT_SPF", 
     "OT_SPF",   "OT_SPF", "OT_CRU",  "OT_MIX_DMF_PEL", "OT_DMF",  "OT_DMF",  "OT_MIX_DMF_PEL",   "OT_DMF",   "OT_DMF",   "OT_DMF",
     "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF",   "SDN_DEM",  "SDN_DEM",  "SDN_DEM",  "SDN_DEM",   "SSC_DEM", "SSC_DEM",  
     "SSC_DEM", "SSC_DEM", "SSC_DEM",  "TBB_CRU",   "TBB_DMF", "TBB_DMF")  

   } else{ 
   if(a_year=="2008"){
     levels(tacsatp$LE_MET) <-   c(  ## REPLACE LEVELS WITH CAUTION ## adapt to your own list!!
     "DRB_MOL", "OT_DMF", "NA", "OT_CRU", "OT_CRU", "OT_CRU", "OT_CRU", "OT_CRU", "OT_CRU",  "OT_MIX_DMF_PEL",    
     "OT_DMF", "OT_DMF", "OT_DMF", "OT_MIX_DMF_PEL", "OT_MIX_DMF_PEL", "OT_DMF", "OT_DMF",  "OT_DMF", "OT_DMF",  "OT_SPF",  
     "OT_SPF", "OT_SPF", "OT_SPF",  "OT_CRU",  "OT_MIX_DMF_PEL", "OT_DMF", "OT_DMF", "OT_DMF", "OT_DMF",  "OT_SPF", 
     "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF", "SDN_DEM", "SDN_DEM", "SDN_DEM","SDN_DEM",  "SSC_DEM", "SSC_DEM",  
     "SSC_DEM", "SSC_DEM",  "TBB_CRU",  "TBB_DMF")
   tacsatp$LE_MET         <- as.character(tacsatp$LE_MET)
   } else{ 
   if(a_year=="2009"){
     levels(tacsatp$LE_MET) <-   c(  ## REPLACE LEVELS WITH CAUTION ## adapt to your own list!!
     "DRB_MOL", "NA", "OT_CRU", "OT_CRU", "OT_CRU", "OT_CRU", "OT_CRU", "OT_CRU", "OT_MIX_DMF_PEL", "OT_DMF", "OT_DMF",
     "OT_DMF", "OT_MIX_DMF_PEL", "OT_MIX_DMF_PEL", "OT_DMF", "OT_DMF", "OTB_DMF", "OT_SPF", "OT_SPF", "OT_SPF", 
     "OT_SPF", "OT_SPF", "OT_CRU", "OT_MIX_DMF_PEL", "OT_DMF", "OT_DMF", "OT_MIX_DMF_PEL", "OT_DMF", "OT_SPF", "OT_SPF", 
     "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF", "SDN_DEM", "SDN_DEM", "SDN_DEM", "SDN_DEM",  "SSC_DEM", "SSC_DEM",  
     "SSC_DEM", "SSC_DEM", "TBB_CRU", "TBB_DMF")  
   tacsatp$LE_MET         <- as.character(tacsatp$LE_MET)
   } else{
   if(a_year=="2010"){
     levels(tacsatp$LE_MET) <-   c(  ## REPLACE LEVELS WITH CAUTION ## adapt to your own list!!
     "DRB_MOL", "NA",  "OT_CRU","OT_CRU", "OT_CRU",   "OT_CRU",  "OT_CRU",  
     "OT_CRU",  "OT_MIX_DMF_PEL",     "OT_DMF", "OT_DMF", "OT_DMF",   "OT_DMF", "OT_MIX_DMF_PEL",  
     "OT_MIX_DMF_PEL",   "OT_DMF",  "OT_DMF",   "OT_DMF",  "OT_DMF",  "OT_SPF",  "OT_SPF",  
     "OT_SPF",  "OT_SPF",   "OT_SPF",   "OT_MIX_DMF_PEL",   "OT_CRU",  "OT_MIX_DMF_PEL",     "OT_DMF",
     "OT_DMF",   "OT_DMF",  "OT_SPF",   "OT_SPF",  "OT_SPF",   "OT_SPF",   "SDN_DEM",
     "SDN_DEM",   "SDN_DEM", "SDN_DEM",  "SSC_DEM",   "SSC_DEM",  "TBB_CRU",   "TBB_DMF")  
   tacsatp$LE_MET         <- as.character(tacsatp$LE_MET)
   } else{
    if(a_year=="2011"){
     levels(tacsatp$LE_MET) <-   c(  ## REPLACE LEVELS WITH CAUTION ## adapt to your own list!!
      "DRB_MOL",      "NA",   "OT_CRU", "OT_CRU",   "OT_CRU",  "OT_MIX_DMF_PEL", "OT_DMF", "OT_MIX_DMF_PEL",  
      "OT_MIX_DMF_PEL",   "OT_DMF",  "OT_MIX_NEP",   "OT_MIX_NEP", "OT_MIX_NEP",   "OT_MIX_NEP",  "OT_SPF",  "OT_SPF",  
      "OT_SPF",   "OT_MIX_DMF_PEL",   "OT_SPF",   "OT_CRU",  "OT_MIX_DMF_PEL",     "OT_DMF", "OT_MIX_DMF_PEL",   "OT_MIX_NEP",  
      "OT_MIX_NEP",  "OT_SPF",   "OT_SPF",  "OT_SPF",   "OT_SPF",   "SDN_DEM", "SDN_DEM",   "SDN_DEM",
      "SDN_DEM",  "SSC_DEM", "SSC_DEM",   "TBB_CRU",   "TBB_DMF",     "TBB_DMF")  
   
    }else{
     if(a_year=="2012") {
       levels(tacsatp$LE_MET) <-   c("DRB_MOL", "NA", "OT_CRU", "OT_CRU", "OT_CRU",  "OT_DMF", "OT_DMF", "OT_MIX_DMF_PEL",
        "OT_MIX_DMF_PEL", "OT_DMF", "OT_MIX_NEP", "OT_MIX_NEP",  "OT_MIX_NEP", "OT_MIX_NEP",
        "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF", "OT_MIX_DMF_PEL", "OT_DMF", "OT_MIX_NEP",
        "OT_MIX_NEP",  "OT_SPF", "OT_SPF", "OT_SPF", "OT_SPF", 
        "SDN_DEM", "SDN_DEM",   "SDN_DEM", "SDN_DEM", "SSC_DEM", "SSC_DEM", "TBB_CRU", "TBB_DMF")
    } else{
      if(a_year=="2013") {
       levels(tacsatp$LE_MET) <-   c(
        "DRB_MOL","NA", "OT_CRU", "OT_CRU",
        "OT_CRU", "OT_DMF" ,"OT_DMF","OT_MIX_DMF_PEL",  
        "OT_MIX_DMF_PEL", "OT_DMF", "OT_MIX_NEP","OT_MIX_NEP","OT_MIX_NEP", "OTB_MIX_NEP",
        "OT_SPF","OT_SPF" , "OT_SPF",   "OT_MIX_DMF_PEL",     "OT_SPF",   "OT_MIX_DMF_PEL",    
        "OT_DMF", "OT_MIX_DMF_PEL",   "OT_MIX_NEP",   "OT_SPF",
        "OT_SPF",   "OT_SPF",  "OT_MIX_DMF_PEL",   "OT_SPF",   "OT_SPF",  
        "SDN_DEM", "SDN_DEM",   "SDN_DEM",
        "SDN_DEM",  "SSC_DEM", "SSC_DEM",   "TBB_CRU",   "TBB_DMF")  
     
    } else{
    stop('adapt the BENTHIS metiers for this year')
    }
   }
  }}}}}}}
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
    storeScheme$years <- a_year  
  }

  tacsatp$year      <- format(tacsatp$SI_DATIM, "%Y")
  require(mixtools)
  activity          <- activityTacsat(tacsatp,units="year",analyse.by="LE_GEAR", storeScheme,
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
  
  #- Gear specific fm parameters (parameters tuned with st=c(2,6))
  #- (should be informed for all below towedGears)
  fm        <- list(TBB=0.13,
                    OTB=0.13, # DNK close to straight line (fm at 0)
                    OTT=0.13,
                    PTB=0.13,
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
        try({
        interpolationcHs    <- interpolateTacsat(tacsatpGearVEREF,
                                                 interval= VMS_ping_rate_in_hour*60,   # THE PING RATE IS COUNTRY-SPECIFIC ##
                                                 margin  = round(VMS_ping_rate_in_hour*60*0.1), # i.e. will make disconnected interpolations if interval out of the 50 70min range
                                                 res     = 100,
                                                 method  = "cHs",
                                                 params  = list(fm=fm[[iGr]],distscale=0,sigline=0,st=c(2,6)),   # rmenber that st not in use....
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
        }, silent=TRUE)
        if(class(interpolationcHs)=="try-error") print('error for this interpolation')
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

#for(a_year in c(2005:2013)) {
#print(a_year)

fls <- dir(file.path(outPath, a_year,"interpolated"))
fls <- fls[grep("tacsatSweptArea_", fls)]

lst <- list(); count <- 0
vid_with_errors <- NA
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
  
  if(any(tacsatIntGearVEREF$SWEPT_AREA_KM2>100, na.rm = TRUE) ) {
    print(paste('check for lat long at 0!! for ', iFile))
    vid_with_errors <- c(vid_with_errors, iFile)
    tacsatIntGearVEREF[!is.na(tacsatIntGearVEREF$SWEPT_AREA_KM2) & tacsatIntGearVEREF$SWEPT_AREA_KM2>100, c("SWEPT_AREA_KM2", "SWEPT_AREA_KM2_LOWER", "SWEPT_AREA_KM2_UPPER")] <- NA
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

#}


#-----------------------------------------------------------------------------
# Create the aggregated swept area dataset
# (TO BE DELIVERED BY EACH PARTNER TO THE WP2 COORDINATOR)
#-----------------------------------------------------------------------------

# once you have your three years ready....
outPath   <- "C:/BENTHIS/outputs/"  # PLEASE ADAPT.

library(vmstools)

#- Combine tacsatSweptArea files from 2010-2012
for(iYr in 2010:2012){
#- read data for this year
load(file.path(outPath, iYr, "tacsatSweptArea.RData"))

#- collate
if(iYr == 2010) tacsatSweptAreaTot <- cbind(tacsatSweptArea,SI_YEAR=iYr)
if(iYr != 2010) tacsatSweptAreaTot <- rbind(tacsatSweptAreaTot,cbind(tacsatSweptArea,SI_YEAR=iYr))
}
tacsatSweptArea <- tacsatSweptAreaTot; rm(tacsatSweptAreaTot)

#- add months and days
tacsatSweptArea$SI_DATE <- as.POSIXct(paste(tacsatSweptArea$SI_DATE,    sep=" "), tz="GMT", format="%d/%m/%Y")
tacsatSweptArea$MONTH   <- format(tacsatSweptArea$SI_DATE, "%m") 
tacsatSweptArea$DAY     <- format(tacsatSweptArea$SI_DATE, "%j")  

#- Get outer ranges of latitude and longitude in dataset
xrange <- range(tacsatSweptArea$SI_LONG,na.rm=T)
yrange <- range(tacsatSweptArea$SI_LATI,na.rm=T)
xrange <- c(floor(xrange[1]),ceiling(xrange[2]))
yrange <- c(floor(yrange[1]),ceiling(yrange[2]))
print(xrange); print(yrange)

#- If xrange and yrange are inappropriate, set your own (rounded) xrange and yrange
if(FALSE){
xrange <- c(-10,20) # DEN
yrange <- c(50,66) # DEN
# xrange <- c(0,9) # NLD
# yrange <- c(51,57) # NLD
xrange[1] <- floor(xrange[1]); xrange[2] <- ceiling(xrange[2])
yrange[1] <- floor(yrange[1]); yrange[2] <- ceiling(yrange[2])
}

#- Set grid
resx <- 1/60 #1 minute
resy <- 1/60 #1 minute
grd <- createGrid(xrange,yrange,resx=1/60,resy=1/60,type="SpatialGrid",exactBorder=T)

#- Grid all tacsatSweptArea data
# Convert all tacsat poins first to SpatialPoints
coords <- SpatialPoints(cbind(SI_LONG=tacsatSweptArea$SI_LONG,SI_LATI=tacsatSweptArea$SI_LATI))
idx <- over(coords,grd)
tacsatSweptArea$grID <- idx

#- Remove records that are not in the study area
tacsatSweptArea <- subset(tacsatSweptArea,is.na(grID)==F)

#-1 Aggregate the results by metier and grid ID (aggregate() can be slow: be patient)
aggTacsatSweptArea <- aggregate(tacsatSweptArea[,c("SWEPT_AREA_KM2",
"SWEPT_AREA_KM2_LOWER",
"SWEPT_AREA_KM2_UPPER")],
by=list(tacsatSweptArea$LE_MET,tacsatSweptArea$grID,tacsatSweptArea$SI_YEAR),sum,na.rm=T)
colnames(aggTacsatSweptArea)[1:3] <- c("LE_MET","grID", "Year")

#- Add midpoint of gridcell to dataset
aggResult <- cbind(aggTacsatSweptArea,CELL_LONG=coordinates(grd)[aggTacsatSweptArea$grID,1],
CELL_LATI=coordinates(grd)[aggTacsatSweptArea$grID,2])
save(aggResult,file=file.path(outPath,"AggregatedSweptArea.RData"))


#-2 Aggregate the results by metier and grid ID and BY MONTH(aggregate() can be slow: be patient)
aggTacsatSweptArea2 <- aggregate(tacsatSweptArea[,c("SWEPT_AREA_KM2",
"SWEPT_AREA_KM2_LOWER", "SWEPT_AREA_KM2_UPPER")],
by=list(tacsatSweptArea$LE_MET,tacsatSweptArea$MONTH,tacsatSweptArea$grID,tacsatSweptArea$SI_YEAR),sum,na.rm=T)
colnames(aggTacsatSweptArea2)[1:4] <- c("LE_MET","MONTH", "grID", "Year")

#- Add midpoint of gridcell to dataset
aggResult <- cbind(aggTacsatSweptArea2,CELL_LONG=coordinates(grd)[aggTacsatSweptArea2$grID,1],
CELL_LATI=coordinates(grd)[aggTacsatSweptArea2$grID,2])
save(aggResult,file=file.path(outPath,"AggregatedSweptAreaMonth.RData"))


#-3 Aggregate the results by metier and grid ID and BY DAY(aggregate() can be slow: be patient)
aggTacsatSweptArea3 <- aggregate(tacsatSweptArea[,c("SWEPT_AREA_KM2",
"SWEPT_AREA_KM2_LOWER", "SWEPT_AREA_KM2_UPPER")],
by=list(tacsatSweptArea$LE_MET, tacsatSweptArea$DAY, tacsatSweptArea$grID,tacsatSweptArea$SI_YEAR),sum,na.rm=T)
colnames(aggTacsatSweptArea3)[1:4] <- c("LE_MET","DAY", "grID", "Year")

#- Add midpoint of gridcell to dataset
aggResult <- cbind(aggTacsatSweptArea3, CELL_LONG=coordinates(grd)[aggTacsatSweptArea3$grID,1],
CELL_LATI=coordinates(grd)[aggTacsatSweptArea3$grID,2])
save(aggResult,file=file.path(outPath,"AggregatedSweptAreaDay.RData"))


#-----------------------------------------------------------------------------
# Create the "missing effort" dataset
# (TO BE DELIVERED BY EACH PARTNER TO THE WP2 COORDINATOR)
#-----------------------------------------------------------------------------

#- Code to get the missing effort in percentages per ICES rectangle
# i.e. the total effort in eflalo compared to the total effort in eflalo from
# VMS-equipped vessels

library(vmstools)

# once you have your three years ready....
dataPath  <- "C:/BENTHIS/EflaloAndTacsat/"


aggResult <- NULL

#- load eflalo 2010-2012
for(iYr in 2010:2012){
  #- read data for this year
  load(file.path(outPath,iYr,"cleanEflalo.RData"))
  load(file.path(outPath,iYr,"tacsatSweptArea.RData"))

  #- Load the datasets
  #data(europa)

  #- Convert time stamps to posixct formats
  eflalo$LE_CDATIM  <- as.POSIXct(eflalo$LE_CDAT,format="%d/%m/%Y",tz="GMT")

  #- Calculate the effort (INTVDAY) in eflalo
  eflalo$INTV       <- c(difftime(eflalo$FT_LDATIM,eflalo$FT_DDATIM,units="mins"))
  eflalo$dummy      <- 1
  eflalo            <- merge(eflalo,aggregate(eflalo$dummy,by=list(eflalo$FT_REF,eflalo$LE_CDATIM),FUN=sum,na.rm=T),by.x=c("FT_REF","LE_CDATIM"),by.y=c("Group.1","Group.2"),all.x=T)
  colnames(eflalo)[length(colnames(eflalo))] <- "NR_FT_REF"
  eflalo$INTVDAY    <- eflalo$INTV / eflalo$NR_FT_REF

  # - we need to retrieve the BENTHIS metiers if not already informed in eflalo.
  # this following object is in the workflow line 347, so partners should have it.
  # if not, then they should apply the same procedure you did to assign BENTHIS metier
  ctry <- "DEN"
  if(ctry=="DEN"){
    load(file=file.path(outPath,iYr,"initVersusBenthisMetiers.RData"))
    eflalo$LE_MET_BENTHIS <- initVersusBenthisMetiers[match(eflalo$LE_MET, initVersusBenthisMetiers$LE_MET_init),'LE_MET']
  }

  # note that the metiers at NAs correspond to all the metiers which are not BENTHIS e.g., gillnets, pots, etc.

  
  # - partner might adapt here if needed (the goal is to subset eflalo for the VMS-equipped vessels only)
  fls          <- dir(file.path(outPath,iYr,"interpolated"))
  vms_equipped <- unique(sapply(fls, function (x) unlist(strsplit(x, "_"))[2]))  
  eflalo_vms   <- eflalo[eflalo$VE_REF %in%  vms_equipped,]

   
  # - do the aggregation
  aggResult_vms<- aggregate(eflalo_vms$INTVDAY, list(eflalo_vms$LE_RECT, eflalo_vms$LE_MET_BENTHIS), sum, na.rm=TRUE)
  colnames(aggResult_vms) <- c('LE_RECT', 'LE_MET', 'INTVDAY')
  aggResult_tot<- aggregate(eflalo$INTVDAY, list(eflalo$LE_RECT,  eflalo$LE_MET_BENTHIS), sum, na.rm=TRUE)
  colnames(aggResult_tot) <- c('LE_RECT', 'LE_MET', 'INTVDAY_TOT')
  aggResult    <- rbind.data.frame(aggResult, cbind.data.frame(merge(aggResult_vms, aggResult_tot), SI_YEAR=iYr))
}




save(aggResult,file=file.path(outPath,paste("missingEffortTable.RData", sep='')))


