`sortTacsat` <-
function(dat){
require(doBy)

if(!"SI_DATIM" %in% colnames(dat)) dat$SI_DATIM  <- as.POSIXct(paste(dat$SI_DATE,  dat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")

  #Sort the tacsat data first by ship, then by date
if("VE_REF" %in% colnames(dat)) dat <- orderBy(~VE_REF+SI_DATIM,data=dat)
if("OB_REF" %in% colnames(dat)) dat <- orderBy(~OB_REF+SI_DATIM,data=dat)

return(dat)}

                                                