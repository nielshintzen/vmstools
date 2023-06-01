#' Find destination from point of origin given bearing and distance
#' 
#' Find destination from point of origin given bearing and distance
#' 
#' 
#' @param lon Longitude of origin
#' @param lat Latitude of origin
#' @param bearing Bearing to destination
#' @param distance Distance to cover to destination in km
#' @return Returnes the destination point(s) as matrix
#' @author Niels T. Hintzen
#' @seealso \code{\link{addWidth}}
#' @references EU Lot 2 project
#' @examples
#' 
#' res <- destFromBearing(rep(2.5,10),rep(51.5,10),
#'                        runif(10,0,360),runif(10,0,0.1))
#' plot(res[,1],res[,2])
#' points(2.5,51.5,cex=1.1,pch=19)
#' 
#' @export destFromBearing
destFromBearing <- function(lon,lat,bearing,distance){

                  pd    <- pi/180
                  dist  <- (distance/6371)
                  y1    <- lat * pd
                  x1    <- lon * pd
                  bear  <- bearing * pd

                  y2  <- asin(sin(y1) * cos(dist) + cos(y1) * sin(dist) * cos(bear))
                  x2  <- x1 + atan2(sin(bear) * sin(dist) * cos(y1),cos(dist) - sin(y1) * sin(y2))

                  x2  <- (x2 + 3*pi) %% (2*pi) - pi
                  y2  <- y2 / pd
                  x2  <- x2 / pd
                return(cbind(x2,y2))}
