ICESarea <- function(tacsat,areas,proj4string=NULL,fast=FALSE){
  require(sp)
  if(!class(areas) %in% c("SpatialPolygons","SpatialPolygonsDataFrame")) stop("'areas' must be specified as class 'SpatialPolygons' or 'SpatialPolygonsDataFrame'")
  if(class(areas) == "SpatialPolygonsDataFrame") areas <- as(areas,"SpatialPolygons")
  #filter NA values
  NAS         <- which(is.na(tacsat$SI_LONG)==F & is.na(tacsat$SI_LATI)==FALSE)
  
  totres      <- rep(NA,length(tacsat$SI_LONG))
  nms         <- unlist(lapply(areas@polygons,function(x){return(x@ID)}))
  
  if(is.null(proj4string)==TRUE & is.na(proj4string(areas))){
    #No projection string used
    spPoint           <- SpatialPoints(data.frame(x=tacsat$SI_LONG[NAS], y=tacsat$SI_LATI[NAS]))
    
    if(nrow(tacsat)>1e4 & fast == FALSE) {
      idx     <- numeric()
      cat(paste("Calculation memory demanding, therefore a stepwise approach is taken. use: fast=TRUE for memory demaning approach\n"))
      for (chunk in 1: ceiling(nrow(tacsat) / 1e4)) # if number of points exceeds 1e4 then chunk the process to avoid lacking of memory from the over() function!
      {
        #cat(paste("chunk...[",(1+(1e4*(chunk-1))),", ", min((1e4*chunk), nrow(tacsat)),"]\n"))
        idx <- c(idx, over(spPoint [(1+(1e4*(chunk-1))): min((1e4*chunk), nrow(tacsat)),] , areas) )
      }
    } else {
      idx               <- over(spPoint, areas)
    }
    
    totres[NAS]       <- nms[idx]
    
  } else {
    #Use projection string
    if(is.null(proj4string))
      proj4string     <- proj4string(areas)
    spPoint           <- SpatialPoints(data.frame(x=tacsat$SI_LONG[NAS],y=tacsat$SI_LATI[NAS]), proj4string=CRS(proj4string))
    if(nrow(tacsat)>1e4 & fast == FALSE) {
      idx     <- numeric()
      cat(paste("Calculation memory demanding, therefore a stepwise approach is taken. use: fast=TRUE for memory demaning approach\n"))
      for (chunk in 1: ceiling(nrow(tacsat) / 1e4)) # if number of points exceeds 1e4 then chunk the process to avoid lacking of memory from the over() function!
      {
        #cat(paste("chunk...[",(1+(1e4*(chunk-1))),", ", min((1e4*chunk), nrow(tacsat)),"]\n"))
        idx <- c(idx, over(spPoint [(1+(1e4*(chunk-1))): min((1e4*chunk), nrow(tacsat)),] , areas) )
      }
    } else {
      idx               <- over(spPoint, areas)
    }
    
    totres[NAS]       <- nms[idx]
  }
  return(totres)}
