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
pointOnLand <- function(tacsat,lands,proj4string=NULL){
      if(class(lands) != "SpatialPolygons") stop("'lands' must be specified as class 'SpatialPolygons'")

      totres      <- rep(0,length(tacsat$SI_LONG))

      if(is.null(proj4string)==TRUE){
        #No projection string used
        if(is.na(proj4string(lands))==FALSE) stop("Projection defined for lands, use proj4string argument in function")
        spPoint           <- SpatialPoints(data.frame(x=tacsat$SI_LONG,y=tacsat$SI_LATI))
        idx               <- over(spPoint,lands)
        totres[which(is.na(idx)==FALSE)] <- 1
        totres[which(is.na(idx)==TRUE)] <- 0

      } else {
          #Use projection string
          proj4string(lands)<- proj4string
          spPoint           <- SpatialPoints(data.frame(x=tacsat$SI_LONG,y=tacsat$SI_LATI),proj4string=proj4string)
          idx               <- over(spPoint,lands)
          totres[which(is.na(idx)==FALSE)] <- 1
          totres[which(is.na(idx)==TRUE)] <- 0
        }
return(totres)}


