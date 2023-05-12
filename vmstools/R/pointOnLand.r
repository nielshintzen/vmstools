#' Find points on land given a set of coordinates
#' 
#' Find the points that are on land given a set of coordinates and polygons
#' that determine the part that is land
#' 
#' With many coordinates, the checking might take longer.
#' 
#' @param tacsat Tacsat file
#' @param lands Polygon of area that is considered to be land
#' @param proj4string Projection string, default to NULL.
#' @return Returns a vector with values 0 and 1. 1 indicating points on land, 0
#' indicating points not on land.
#' @author Niels T. Hintzen
#' @seealso \code{\link{pointInHarbour}}
#' @references EU Lot 2 project
#' @examples
#' 
#' 
#' data(tacsat)
#' data(europa)
#' tacsat  <- tacsat[1:1000,]
#' tacsat  <- sortTacsat(tacsat)
#' 
#' pols    <- lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),
#'                 function(x){data.frame(SI_LONG=subset(europa,SID==x)$X,
#'                                        SI_LATI=subset(europa,SID==x)$Y)}))
#' idx     <- pointOnLand(tacsat,pols);
#' idx     <- which(idx == 1)
#' 
#' plotMap(europa,xlim=c(0,10),ylim=c(48,62))
#' points(tacsat$SI_LONG[idx],tacsat$SI_LATI[idx],col="red",cex=0.5,pch=19)
#' 
#' @export pointOnLand
pointOnLand <- function(tacsat,lands,st_crs=NULL){
      if(!"sf" %in% class(lands)) stop("'lands' must be specified as class 'sf'")

      totres      <- rep(0,length(tacsat$SI_LONG))

      if(is.null(st_crs)==TRUE){
        #No projection string used
        if(is.na(st_crs(lands))==FALSE) stop("Projection defined for lands, use st_crs argument in function")
        spPoint           <- st_as_sf(tacsat,coords=c("SI_LONG","SI_LATI"))
        idx               <- st_over(spPoint,lands)
        totres[which(is.na(idx)==FALSE)] <- 1
        totres[which(is.na(idx)==TRUE)] <- 0

      } else {
          #Use projection string
          if(is.na(st_crs(lands))){
            st_crs(lands)     <- st_crs
          } else {
            lands           <- st_transform(lands,crs=st_crs)
          }

          spPoint           <- st_as_sf(tacsat,coords=c("SI_LONG","SI_LATI"),crs=st_crs)
          idx               <- st_over(spPoint,lands)
          totres[which(is.na(idx)==FALSE)] <- 1
          totres[which(is.na(idx)==TRUE)] <- 0
        }
return(totres)}
