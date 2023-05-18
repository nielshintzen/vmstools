#' Creates a 'SF polygon' object from longitude-latitude dataset
#' 
#' Takes one set of longitude-lattitude or a list of lon-lat objects and
#' converts it into one sfc object.
#' 
#' 
#' @param SI_LONG Set with longitude values
#' @param SI_LATI Set with latitude values
#' @param lst list of objects, each with longitude and latitude values
#' @return Returns a 'sfc' object from the package 'sf'
#' @author Niels T. Hintzen
#' @seealso \code{\link{st_polygon}}, \code{\link{st_multipolygon}}
#' @references EU Lot 2 project
#' @examples
#' 
#' data(europa)
#'
#' eurPols     <- lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),
#'                       function(x){data.frame(
#'                          SI_LONG=subset(europa,SID==x)$X,
#'                          SI_LATI=subset(europa,SID==x)$Y)}))
#'
#' area        <- lonLat2SpatialPolygons(SI_LONG=c(2,2.5,2.7,2.1),SI_LATI=c(54,54.2,55.8,55.6))
#' plot(eurPols,col="green",xlim=c(-4,10),ylim=c(48,62)); axis(1);axis(2);box()
#' plot(area,add=TRUE,col="red")
#' 
#' 
#' @export ICESrectangle2CSquare


lonLat2SFPolygons <- function(SI_LONG=NULL,SI_LATI=NULL,lst=NULL){
  if((is.null(SI_LONG)==TRUE | is.null(SI_LATI)==TRUE) & is.null(lst)==TRUE) stop("Specify either longitude and latitude vectors or a list containing a dataframe with specified longitude and latitude vectors")

  #- Return spatial polygons object when SI_LONG and SI_LATI vectors are supplied
  if(is.null(SI_LONG)==FALSE | is.null(SI_LATI)==FALSE){
    if(length(SI_LONG) != length(SI_LATI)) stop("Length SI_LONG not equal to length SI_LATI")
    brks      <- which(is.na(SI_LONG)==TRUE | is.na(SI_LATI) == TRUE)
    storePols <- list()
    #- If no breaks
    if(length(brks)==0){
      if((SI_LONG[1] != SI_LONG[length(SI_LONG)]) | (SI_LATI[1] != SI_LATI[length(SI_LATI)])){
        coords <- st_multipoint(rbind(cbind(SI_LONG,SI_LATI),
                                      cbind(SI_LONG[1],SI_LATI[1])))
      } else {
          coords <- st_multipoint(cbind(SI_LONG,SI_LATI))
        }
      pols         <- st_polygon(list(coords))
      storePols[[1]] <- pols
    } else {
        #- If breaks
        for(i in 1:length(brks)){
          if(i == 1){
            if((SI_LONG[1] != SI_LONG[(brks[i]-1)]) | (SI_LATI[1] != SI_LATI[(brks[i]-1)])){
              coords <- st_multipoint(rbind(cbind(SI_LONG[1:(brks[i]-1)],SI_LATI[1:(brks[i]-1)]),
                                          cbind(SI_LONG[1],SI_LATI[1])))
            } else {
                coords <- st_multipoint(cbind(SI_LONG[1:(brks[i]-1)],SI_LATI[1:(brks[i]-1)]))
              }
          } else {
              if((SI_LONG[(brks[i-1]+1)] != SI_LONG[(brks[i]-1)]) | (SI_LATI[(brks[i-1]+1)] != SI_LATI[(brks[i]-1)])){
                coords <- st_multipoint(rbind(cbind(SI_LONG[(brks[i-1]+1):(brks[i]-1)],SI_LATI[(brks[i-1]+1):(brks[i]-1)]),
                                            cbind(SI_LONG[(brks[i-1]+1)],SI_LATI[(brks[i-1]+1)])))
              } else {
                  coords <- st_multipoint(cbind(SI_LONG[(brks[i-1]+1):(brks[i]-1)],SI_LATI[(brks[i-1]+1):(brks[i]-1)]))
                }
            }
          pols <- st_polygon(list(coords))
          storePols[[i]] <- pols
        }
      }
    spol        <- st_multipolygon(storePols)
  }
  #- Return spatial polygons object when list with SI_LONG,SI_LATI dataframe is given
  if(is.null(lst)==FALSE){
    if(is.null(names(lst))==TRUE) names(lst) <- 1:length(lst)
    storePols <- list()
    counter   <- 1
    for(j in 1:length(lst)){
      brks      <- which(is.na(lst[[j]]$SI_LONG)==TRUE | is.na(lst[[j]]$SI_LATI) == TRUE)
      #- If no breaks
      if(length(brks)==0){
        if((lst[[j]]$SI_LONG[1] != lst[[j]]$SI_LONG[nrow(lst[[j]])]) | (lst[[j]]$SI_LATI[1] != lst[[j]]$SI_LATI[nrow(lst[[j]])])){
          coords <- st_multipoint(rbind(cbind(lst[[j]]$SI_LONG,lst[[j]]$SI_LATI),
                                      cbind(lst[[j]]$SI_LONG[1],lst[[j]]$SI_LATI[1])))
        } else {
            coords <- st_multipoint(cbind(lst[[j]]$SI_LONG,lst[[j]]$SI_LATI))
          }
        pols <- st_polygon(list(coords))
        storePols[[counter]] <- pols
        counter <- counter+1
      } else {
          #- If breaks
          if(length(which(diff(brks)==1)>0)) stop("Only 1 row of NA is allowed in between values of SI_LONG and SI_LATI")
          for(i in 1:length(brks)){
            if(i == 1){
              if((lst[[j]]$SI_LONG[1] != lst[[j]]$SI_LONG[(brks[i]-1)]) | (lst[[j]]$SI_LATI[1] != lst[[j]]$SI_LATI[(brks[i]-1)])){
                coords <- st_multipoint(rbind(cbind(lst[[j]]$SI_LONG[1:(brks[i]-1)],lst[[j]]$SI_LATI[1:(brks[i]-1)]),
                                            cbind(lst[[j]]$SI_LONG[1],lst[[j]]$SI_LATI[1])))
              } else {
                  coords <- st_multipoint(cbind(lst[[j]]$SI_LONG[1:(brks[i]-1)],lst[[j]]$SI_LATI[1:(brks[i]-1)]))
                }
            } else {
                if((lst[[j]]$SI_LONG[(brks[i-1]+1)] != lst[[j]]$SI_LONG[(brks[i]-1)]) | (lst[[j]]$SI_LATI[(brks[i-1]+1)] != lst[[j]]$SI_LATI[(brks[i]-1)])){
                  coords <- st_multipoint(rbind(cbind(lst[[j]]$SI_LONG[(brks[i-1]+1):(brks[i]-1)],lst[[j]]$SI_LATI[(brks[i-1]+1):(brks[i]-1)]),
                                              cbind(lst[[j]]$SI_LONG[(brks[i-1]+1)],lst[[j]]$SI_LATI[(brks[i-1]+1)])))
              } else {
                    coords <- st_multipoint(cbind(lst[[j]]$SI_LONG[(brks[i-1]+1):(brks[i]-1)],lst[[j]]$SI_LATI[(brks[i-1]+1):(brks[i]-1)]))
                }
              }
            pols <- st_polygon(list(coords))
            storePols[[counter]] <- pols
            counter     <- counter + 1
          }
        }
    }
    #- Create spatial polygons
    spol        <- storePols
  }
  return(st_as_sfc(spol))}
        
