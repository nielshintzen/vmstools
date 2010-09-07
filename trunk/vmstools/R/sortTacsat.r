`sortTacsat` <-
function(tacsat){

  #Load a dataframe sorter library      
library(doBy)
SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S")
tacsat    <- cbind(tacsat,SI_DATIM)
  #Sort the tacsat data first by ship, then by date
tacsat <- orderBy(~VE_REF+SI_DATIM,data=tacsat)
return(tacsat[,-grep("SI_DATIM",colnames(tacsat))])}

