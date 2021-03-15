#' Compute distance between two points on a sphere (approximation of the earth)
#' 
#' Compute the distance between two GPS locations defined in longitude and
#' latitude notation on the earth. The earth is assumed to have a perfect
#' spherical shape. Distance is returned in km.
#' 
#' 
#' @param lon Longitude of point 2
#' @param lat Latitude of point 2
#' @param lonRef Longitude of point 1
#' @param latRef Latitude of point 1
#' @author Niels T. Hintzen
#' @seealso \code{\link{km2Degree}}, \code{\link{degree2Km}},
#' \code{\link{lonLatRatio}}
#' @references EU Lot 2 project, based on the Haversine formula, see also:
#' Hintzen et al. 2010 Fisheries Research
#' @examples
#' 
#' lon <- -4
#' lat <- 50
#' lonRef <- -4.2
#' latRef <- 51
#' 
#' distance(lon,lat,lonRef,latRef) #112.09
#' 
#' @export distance
`distance` <-
function(lon,lat,lonRef,latRef){

                    pd <- pi/180

                    a1<- sin(((latRef-lat)*pd)/2)
                    a2<- cos(lat*pd)
                    a3<- cos(latRef*pd)
                    a4<- sin(((lonRef-lon)*pd)/2)
                    a <- a1*a1+a2*a3*a4*a4

                                      c <- 2*atan2(sqrt(a),sqrt(1-a));
                    return(6371*c)}
