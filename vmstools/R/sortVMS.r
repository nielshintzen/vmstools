`sortVMS` <-
function(VMS){

  #Load a dataframe sorter library      
library(doBy)

  #Sort the VMS data first by ship, then by date and alternatively also by speed and heading
VMS <- orderBy(~ship+date+speed+heading,data=VMS)
return(VMS)}

