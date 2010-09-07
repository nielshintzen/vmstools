`sortTacsat` <-
function(tacsat){

  #Load a dataframe sorter library      
library(doBy)
SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S")
  #Sort the tacsat data first by ship, then by date
tacsat <- orderBy(~VI_REF+SI_DATIM,data=tacsat)
return(tacsat)}

