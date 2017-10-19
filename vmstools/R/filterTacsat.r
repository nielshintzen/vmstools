`filterTacsat` <-
function(tacsat
                         ,st=c(2,6)                #Speed threshold points (two values), NULL means use all points
                         ,hd=NULL                  #Heading threshold points (two values), NULL means use all points
                         ,remDup=TRUE              #Specify if you want to remove duplicated VMS records (these should not occur in the first place)
                         ){
  VMS <- tacsat
if(!"SI_DATIM" %in% colnames(VMS))  VMS$SI_DATIM     <- as.POSIXct(paste(VMS$SI_DATE,  VMS$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
  #Remove duplicate records
if(is.null(remDup)==FALSE){
  uniqueVMS <- which(duplicated(VMS[,c("VE_REF","SI_DATIM")])==FALSE)
  VMS.      <- VMS[uniqueVMS,]
  if(dim(VMS.)[1] != dim(VMS)[1]) warning("duplicate records have been removed")
} else {
    VMS. <- VMS
  }

if(is.null(hd)==FALSE){
  warning("It is assumed that VMS data is sorted by vessel and date")
    #Calculate the difference in heading between succeeding VMS datapoints
  diffHeading                           <- abs(VMS.$SI_HE[2:dim(VMS.)[1]] - VMS.$SI_HE[1:(dim(VMS.)[1]-1)])
    #If there is a change to a new vessel in the data, put diffHeading to NA
  diffHeading[unlist(lapply(as.list(unique(VMS.$VE_REF)),function(x){return(which(VMS.$VE_REF == x)[1])}))[-1]-1] <- NA
  VMS.$diffHeading                      <- 0
    #Store the difference in heading
  VMS.$diffHeading[1:(dim(VMS.)[1]-1)]  <- diffHeading
    #Subset based on the absolute difference in heading
  VMS.                                  <- subset(VMS.,diffHeading %in% seq(hd[1],hd[2],1))
}

if(is.null(st)==FALSE){
    #Subset based on the speed ranges
  if(length(st)>2) stop("More than two speed values selected")
  VMS. <- VMS.[which(VMS.$SI_SP >= st[1] & VMS.$SI_SP <= st[2]),]
}
return(VMS.[,-grep("SI_DATIM",colnames(VMS.))])}

