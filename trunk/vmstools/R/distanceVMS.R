distanceVMS <- function(VMS,index){

                            VMS[index[1,1]:index[1,2],]

                   res <- unlist(lapply(as.list(1:dim(index)[1]),function(x){
                              iS  <- index[x,1]
                              iE  <- index[x,2]
                              iL  <- iE-iS+1
                              res <- distance(VMS[iS:iE,]$declon[2:iL],VMS[iS:iE,]$declat[2:iL],VMS[iS:iE,]$declon[1:(iL-1)],VMS[iS:iE,]$declat[1:(iL-1)])
                           return(sum(res,na.rm=T))}))

               return(res)}