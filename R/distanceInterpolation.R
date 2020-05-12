distanceInterpolation <- function(interpolation){
                                               
                            res <- unlist(lapply(interpolation,function(x){
                                                 dims        <- dim(x)
                                                 res         <- distance(x[3:dims[1],1],x[3:dims[1],2],x[2:(dims[1]-1),1],x[2:(dims[1]-1),2])
                                                 return(sum(res,na.rm=TRUE))}))
                     
                            return(res)}


  


  