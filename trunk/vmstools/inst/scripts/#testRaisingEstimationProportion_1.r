library(vmstools)
data(tacsat)
data(eflalo)

tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp <- sortTacsat(tacsatp)
tacsatp <- intervalTacsat(tacsatp,level="trip",fill.na=T,weight=c(0.5,0.5))
tacsatp1<- subset(tacsatp,FT_REF==0)
tacsatp2<- subset(tacsatp,FT_REF!=0)
tacsatp1<- intervalTacsat(tacsatp1,level="vessel",fill.na=T,weight=c(0.5,0.5))
tacsatp <- rbind(tacsatp1,tacsatp2)
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]
tacsatp <- filterTacsat(tacsatp,st=c(4,8))
tacsatp$INTV[which(tacsatp$INTV>232)] <- 115
tacsatp$INTV[which(is.na(tacsatp$INTV)==T)] <- 115
tacsat  <- tacsatp
tacsat$SI_STATE <- 1

source("D:/Repository/VMStools/vmstools/R/estimatePropFishing.r")
source("D:/Repository/VMStools/vmstools/R/raiseTacsat.r")

eflaloOrig <- eflalo
tacsatOrig <- tacsat
eflalo      <- eflaloOrig
tacsat      <- tacsatOrig

eflalo  <- estimatePropFishing(tacsat,eflalo,by=c("VE_REF","FT_REF"))
tacsat  <- raiseTacsat(tacsat,eflalo,by=c("LE_GEAR","SI_DAY","SI_YEAR"),sortBy=F)

