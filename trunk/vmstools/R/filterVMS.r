`filterVMS` <-
function(VMS
                         ,st=c(2,6)                #Speed threshold points (two values), NULL means use all points
                         ,hd=NULL                  #Heading threshold points (two values), NULL means use all points
                         ,remDup=T                 #Specify if you want to remove duplicated VMS records (these should not occur in the first place)
                         ){                            
  #Remove duplicate records
if(is.null(remDup)==F){
  uniqueVMS <- which(duplicated(VMS[,c("ship","date")])==F)
  VMS.      <- VMS[uniqueVMS,]
  if(dim(VMS.)[1] != dim(VMS)[1]) warning("duplicate records have been removed")
} else {
    VMS. <- VMS
  }

if(is.null(hd)==F){
  warning("It is assumed that VMS data is sorted by vessel and date")
    #Calculate the difference in heading between succeeding VMS datapoints
  diffHeading                           <- abs(VMS.$heading[2:dim(VMS)[1]] - VMS.$heading[1:(dim(VMS)[1]-1)])
    #If there is a change to a new vessel in the data, put diffHeading to NA
  diffHeading[unlist(lapply(as.list(unique(VMS.$ship)),function(x){return(which(VMS.$ship == x)[1])}))[-1]-1] <- NA
  VMS.$diffHeading                      <- 0
    #Store the difference in heading
  VMS.$diffHeading[1:(dim(VMS.)[1]-1)]  <- diffHeading
    #Subset based on the absolute difference in heading
  VMS.                                  <- subset(VMS.,diffHeading %in% seq(hd[1],hd[2],1))
}

if(is.null(st)==F){
    #Subset based on the speed ranges
  if(length(st)>2) stop("More than two speed values selected")
  VMS. <- subset(VMS.,speed %in% seq(st[1],st[2],0.1))
}
return(VMS.)}

