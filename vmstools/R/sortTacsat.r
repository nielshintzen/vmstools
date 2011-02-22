`sortTacsat` <-
function(tacsat){

if(!"SI_DATIM" %in% colnames(tacsat)) tacsat$SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")

  #Sort the tacsat data first by ship, then by date
tacsat <- orderBy(~VE_REF+SI_DATIM,data=tacsat)

return(tacsat)}

                                                