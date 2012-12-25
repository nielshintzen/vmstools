distanceTacsat <- function(tacsat,index){

                   res <- unlist(lapply(as.list(1:dim(index)[1]),function(x){
                              iS  <- index[x,1]
                              iE  <- index[x,2]
                              iL  <- iE-iS+1
                              res <- distance(tacsat[iS:iE,]$SI_LONG[2:iL],tacsat[iS:iE,]$SI_LATI[2:iL],
                              tacsat[iS:iE,]$SI_LONG[1:(iL-1)],tacsat[iS:iE,]$SI_LATI[1:(iL-1)])
                           return(sum(res,na.rm=TRUE))}))

               return(res)}