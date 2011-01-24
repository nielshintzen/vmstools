\name{VMStools-package}
\alias{VMStools-package}
\alias{VMStools}
\docType{package}
\title{
development of R scripts to setup the workflow of vms and logbook data to compute pressure indicators
}
\description{
-Design of routines to import VMS and Logbook data\\*
-Design of routines to cluster metiers\\*
-Design of routines to link VMS and Logbook data\\*
-Design of routines to define vessel activity \\*
-Design of routines to interpolate VMS data \\*
-Design of routines to compute pressure indicators
}
\details{
\tabular{ll}{
Package: \tab VMStools\cr
Type: \tab Package\cr
Version: \tab 0.42\cr
Date: \tab 2011-01-24\cr
License: \tab \cr
LazyLoad: \tab yes\cr
}
data(tacsat)
data(eflalo)
data(harbours)

#-------------------------------------------------------------------------------
#- 2) Clean the tacsat data
#-------------------------------------------------------------------------------
  #- Remove points that cannot be possible
idx             <- which(abs(tacsat$SI_LATI) > 90 | abs(tacsat$SI_LONG) > 180)
idx             <- unique(c(idx,which(tacsat$SI_HE < 0 | tacsat$SI_HE > 360)))
idx             <- unique(c(idx,which(tacsat$SI_SP > spThres)))
tacsat          <- tacsat[-idx,]

  #- Remove duplicate records
tacsat$SI_DATIM <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S")
uniqueTacsat    <- paste(tacsat$VE_REF,tacsat$SI_LATI,tacsat$SI_LONG,tacsat$SI_DATIM)
tacsat          <- tacsat[!duplicated(uniqueTacsat),]

  #- Remove points in harbour
idx             <- pointInHarbour(tacsat$SI_LONG,tacsat$SI_LATI,harbours); pih <- tacsat[which(idx == 1),]
save(pih,file=paste(outPath,"pointInHarbour.RData",sep=""))
tacsat          <- tacsat[which(idx == 0),]

  #- Remove points on land
europe          <- list()
for(iLand in c("Netherlands","Germany","Belgium","France","UK","Denmark","Ireland","Norway","Sweden","Iceland")){
  europe[[iLand]] <- map("worldHires",iLand,main=iLand,plot=F,res=0); #title(iLand)
}
europe[["Norway"]]$y[which(europe[["Norway"]]$y < 50)]  <- NA
europe[["UK"]]$y    [which(europe[["UK"]]$y < 45)]      <- NA; europe[["UK"]]$x    [which(europe[["UK"]]$x < -10)]          <- NA
europe[["France"]]$y[which(europe[["France"]]$y < 40)]  <- NA; europe[["France"]]$x[which(abs(europe[["France"]]$x) > 10)]  <- NA

idx             <- pointOnLand(tacsat,europe); pol <- tacsat[which(idx == 1),]
tacsat          <- tacsat[which(idx == 0),]

  #- Remove points which are pseudo duplicates as they have an interval rate < x minutes
tacsat          <- sortTacsat(tacsat)
tacsatp         <- intervalTacsat(tacsat,level="vessel")
tacsat          <- tacsatp[which(tacsatp$INTV > intThres),-grep("INTV",colnames(tacsatp))]

#-------------------------------------------------------------------------------
#- 3) Clean the eflalo data
#-------------------------------------------------------------------------------
  #- Remove eflalo data with mesh size that are too high
idx             <- which(eflalo$LE_MSZ >= msThres)
eflalo          <- eflalo[-idx,]

  #- Remove outlying catch records
specs           <- substr(colnames(eflalo)[grep("KG",colnames(eflalo))],7,9)
  #- Calculate the landing boundaries
specBounds      <- lapply(as.list(specs),function(x){ idx   <- grep(x,colnames(eflalo))[grep("KG",colnames(eflalo)[grep(x,colnames(eflalo))])];
                                                        wgh   <- sort(unique(eflalo[which(eflalo[,idx]>0),idx]));
                                                        difw  <- diff(log10(wgh));
                                                        return(ifelse(any(difw > lanThres),wgh[rev(which(difw <= lanThres))],ifelse(length(wgh)==0,0,max(wgh,na.rm=T))))})
specBounds      <- cbind(specs,unlist(specBounds))
  #- Get the index of each of the species
idx             <- unlist(lapply(as.list(specs),function(x){ idx   <- grep(x,colnames(eflalo))[grep("KG",colnames(eflalo)[grep(x,colnames(eflalo))])];
                                                        return(idx)}))
  #- If landing > bound -> NA
for(iSpec in idx) eflalo[which(eflalo[,iSpec] > an(specBounds[(iSpec-idx[1]+1),2])),iSpec] <- NA

  #- Remove non-unique trip numbers
eflalo <- eflalo[!duplicated(eflalo$LE_ID),]

  #- Remove records with arrival date before departure date
eflalop           <- eflalo
eflalop$FT_DDATIM <- as.POSIXct(paste(eflalo$FT_DDAT,  eflalo$FT_DTIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S")
eflalop$FT_LDATIM <- as.POSIXct(paste(eflalo$FT_LDAT,  eflalo$FT_LTIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S")
idx               <- which(eflalop$FT_LDATIM >= eflalop$FT_DDATIM)
eflalo            <- eflalo[idx,]

#-------------------------------------------------------------------------------
#- 4) Merge the tacsat and eflalo data together
#-------------------------------------------------------------------------------

tacsatp           <- mergeEflalo2Tacsat(eflalo,tacsat)

#-------------------------------------------------------------------------------
#- 5) Calculate effort in both tacsat as eflalo
#-------------------------------------------------------------------------------

teffort           <- effort(tacsatp,level="trip",unit="hours")
eeffort           <- effort(eflalo,level="trip",unit="hours")

#-------------------------------------------------------------------------------
#- 6) Define activitity
#-------------------------------------------------------------------------------

uniqueTypes   <- unique(eflalo[,c("FT_REF","LE_GEAR")])
Fock2008      <- data.frame(gear=   c("DRB","FPO","GN","GNS","LHP","OTM","PTM","OTB","PTB","SDN","SPR","SSC","SX","TBS","TBB","TX","NON"),
                            speeds= c(2.1,2.8,1.7,2,1.32,2.8,3.3,3,3.4,1.28,2.5,1.51,2.1,3.2,6,3.2,4.2))
replaceTypes  <- data.frame(gear=   c("GND","BNT","LL","LLS","LLD","LHM","BTF","GTR","LH","HMD","OTG","GTN","PS","MIS","OTT","NVT","NK","OFG","AG","TGB"),
                            replace=c("GN","GN","LHP","LHP","LHP","LHP","LHP","GN","LHP","DRB","TX","GN","NON","NON","OTM","NON","NON","NON","NON","NON"))

  #- Replace the gears that are not in Fock2008 with ones that are represented by Fock2008 according to the replaceTypes table
for(iGear in 1:nrow(replaceTypes)) eflalo$LE_GEAR <- replace(eflalo$LE_GEAR,which(eflalo$LE_GEAR == ac(replaceTypes[iGear,"gear"])),ac(replaceTypes[iGear,"replace"]))

  #- Replace TBB gears with hp < 300 to TBS
eflalo$LE_GEAR[which(eflalo$LE_GEAR == "TBB" & (eflalo$VE_KW/1.34102209) < 300)] <- "TBS"

  #- Merge gears from eflalo to tacsat. Note that in some trips more than 1 gear type is used. Here only the first one counts
tacsatp$LE_GEAR   <- merge(tacsatp[,c("VE_REF","FT_REF")],uniqueTypes[which(!duplicated(uniqueTypes$FT_REF)),],by="FT_REF",all.x=T)$LE_GEAR

  #- Define activity
tacsatp$SI_STATE  <- 0
for(iGear in 1:nrow(Fock2008)){
  idx         <- which(tacsatp$LE_GEAR == Fock2008$gear[iGear])
  if(length(idx)>0) idx                   <- idx[which(tacsatp$SI_SP[idx] >= ((1-relSpThres)*Fock2008$speeds[iGear]) & tacsatp$SI_SP[idx] <= ((1+relSpThres)*Fock2008$speeds[iGear]))]
  if(length(idx)>0) tacsatp$SI_STATE[idx] <- 1
}

#-------------------------------------------------------------------------------
#- 7) Dispatch landings at the ping scale
#-------------------------------------------------------------------------------

  #- Choose species
  specs <- c("HER","PLE","SOL")
  idxkg <- grep("KG",colnames(eflalo))
  idxeu <- grep("EURO",colnames(eflalo))
  idx   <- c(idxkg[which(!substr(colnames(eflalo)[idxkg],7,9) %in% specs)],idxeu[which(!substr(colnames(eflalo)[idxeu],9,11) %in% specs)])
tacsatEflalo <- splitAmongPings(tacsat=tacsatp,eflalo=eflalo[,-idx],variable="all",level="day",conserve=T)

#-------------------------------------------------------------------------------
#- 8) Interpolate the tacsat data + dispatch at the ping scale
#-------------------------------------------------------------------------------

  #- Prepare tacsat file for interpolation by removing all points without state == 1
tacsatInter         <- subset(tacsatp,SI_STATE == 1 & FT_REF != 0)
tacsatInter         <- sortTacsat(tacsatInter)

tacsatInterShip     <- split(tacsatInter,af(tacsatInter$VE_REF))
interpolation       <- list()
for(iShip in 1:length(tacsatInterShip))print(iShip); interpolation[[iShip]] <- interpolateTacsat(tacsatInterShip[[iShip]],interval=120,margin=12,res=100,method="cHs",params=list(fm=0.45,distscale=20,sigline=0.2,st=c(4,8)),headingAdjustment=0)
tacsatInterList     <- list()
for(iShip in 1:length(interpolation)){
  if(class(try(nrow(interpolation[[iShip]][[1]])))!="try-error"){
    print(iShip)
    tacsatInterList[[iShip]]   <- interpolation2Tacsat(interpolation[[iShip]],tacsatInterShip[[iShip]],npoints=10)
  }
}
tacsatInter         <- do.call(rbind,tacsatInterList)
tacsatInter$SI_STATE<- 1
tacsatInter         <- sortTacsat(tacsatInter)
tacsatInterEflalo   <- splitAmongPings(tacsat=tacsatInter,eflalo=eflalo[,-idx],variable="all",level="day",conserve=T)

#-------------------------------------------------------------------------------
#- 9) Calculate indicators
#-------------------------------------------------------------------------------
tacsatp             <- sortTacsat(tacsatp)
tacsatp$SI_INTV     <- intervalTacsat(tacsatp,level="vessel")$INTV
tacsatp$SI_INTV[which(tacsatp$SI_INTV > 240)] <- NA
ind5                <- DCFIndicator5(tacsatp,minThreshold=10,cellresX=0.1,cellresY=0.05,calcAreaMethod="Trapezoid",plotMapTF=FALSE,exportTableName="")
ind6                <- DCFIndicator6(tacsatp,pctThreshold=90,cellresX=0.1,cellresY=0.05,calcAreaMethod="Trapezoid",plotMapTF=FALSE,exportTableName="")

}
\author{
Niels T. Hintzen, Doug Beare, Francois Bastardie, Nicolas Deporte, Andy South, Neil Campbell, Fabrizio Manco, Clara Ulrich Rescan,
Stephanie Mahevas, Sebastien Demaneche, Josefine Egekvist, Hans Gerritsen, Matthew Parker-Humphreys

Maintainer: Niels T. Hintzen <niels.hintzen@wur.nl>, Francois Basterdie <fba@aqua.dtu.dk>
}
\references{EU Lot 2 project
}
