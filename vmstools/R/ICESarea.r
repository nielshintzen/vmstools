#' Get ICES area from coordinates
#' 
#' Get the ICES area from any lon,lat position, given that this position is
#' within the ICES region.
#' 
#' 
#' @param tacsat dataframe given that they have 'SI_LONG' and 'SI_LATI' columns
#' (either tacsat format or other dataset with SI_LONG and SI_LATI columns)
#' @param areas ICES areas as SpatialPolygons
#' @param proj4string Projection string, default to NULL.
#' @param fast If memory allocation is not a problem, a faster version can be
#' switched on
#' @return Returns the areas as a vector
#' @author Niels T. Hintzen
#' @seealso \code{\link{ICESrectangle}}, \code{\link{ICESrectangle2LonLat}}
#' @references EU Lot 2 project
#' @examples
#' 
#' data(ICESareas)
#' res   <- data.frame(SI_LONG = c(1,2,2,4,2),
#'                     SI_LATI = c(53,53.2,54,56.7,55.2))
#' areas <- ICESarea(res,ICESareas,st_crs=4326)
#' ICESareas$Area_Full[areas]
#' 
#' @export ICESarea
ICESarea <- function(tacsat,areas,st_crs=NULL,fast=FALSE){
  require(sf)
  if(!"sf" %in% class(areas)) stop("'areas' must be specified as class 'sf'")
  #filter NA values
  NAS         <- which(is.na(tacsat$SI_LONG)==F & is.na(tacsat$SI_LATI)==FALSE)
  
  totres      <- rep(NA,length(tacsat$SI_LONG))
  nms         <- st_drop_geometry(areas)$OBJECTID
  
  if(is.null(st_crs)==TRUE & is.na(st_crs(areas))){
    #No projection string used
    spPoint           <- st_as_sf(tacsat[NAS,],coords=c("SI_LONG","SI_LATI"))
    
    if(nrow(tacsat)>1e4 & fast == FALSE) {
      idx     <- numeric()
      cat(paste("Calculation memory demanding, therefore a stepwise approach is taken. use: fast=TRUE for memory demaning approach\n"))
      for (chunk in 1: ceiling(nrow(tacsat) / 1e4)) # if number of points exceeds 1e4 then chunk the process to avoid lacking of memory from the over() function!
      {
        #cat(paste("chunk...[",(1+(1e4*(chunk-1))),", ", min((1e4*chunk), nrow(tacsat)),"]\n"))
        idx <- c(idx,sapply(st_intersects(spPoint [(1+(1e4*(chunk-1))): min((1e4*chunk), nrow(tacsat)),] ,areas), function(z) if (length(z)==0) NA_integer_ else z[1]))
      }
    } else {
      idx               <- sapply(st_intersects(spPoint, areas), function(z) if (length(z)==0) NA_integer_ else z[1])
    }
    totres[NAS]       <- nms[idx]
    
  } else {
    #Use projection string
    if(is.null(st_crs))
      st_crs     <- st_crs(areas)
    spPoint      <- st_as_sf(tacsat[NAS,],coords=c("SI_LONG","SI_LATI"),crs=st_crs)

    if(nrow(tacsat)>1e4 & fast == FALSE) {
      idx     <- numeric()
      cat(paste("Calculation memory demanding, therefore a stepwise approach is taken. use: fast=TRUE for memory demaning approach\n"))
      for (chunk in 1: ceiling(nrow(tacsat) / 1e4)) # if number of points exceeds 1e4 then chunk the process to avoid lacking of memory from the over() function!
      {
        #cat(paste("chunk...[",(1+(1e4*(chunk-1))),", ", min((1e4*chunk), nrow(tacsat)),"]\n"))
        idx <- c(idx,sapply(st_intersects(spPoint [(1+(1e4*(chunk-1))): min((1e4*chunk), nrow(tacsat)),] ,areas), function(z) if (length(z)==0) NA_integer_ else z[1]))
      }
    } else {
      idx               <- sapply(st_intersects(spPoint, st_geometry(areas)), function(z) if (length(z)==0) NA_integer_ else z[1])
    }
    
    totres[NAS]       <- nms[idx]
  }
  return(totres)}

