#' Compute distance from degrees into kilometers
#' 
#' Function transformes the distance expressed in degrees into kilometers. This
#' based on the GPS location of a point.
#' 
#' 
#' @param lon Longitude of the GPS position
#' @param lat Latitude of the GPS positiona
#' @param degree Value in degrees to turn into Km
#' @note Computation of Km is approximation based on the Haversine formula
#' @author Niels T. Hintzen
#' @seealso \code{\link{distance}}, \code{\link{km2Degree}},
#' \code{\link{lonLatRatio}}
#' @references EU lot 2 project
#' @examples
#' 
#' lon <- -4
#' lat <- 50
#' degree <- 1.601833
#' 
#' degree2Km(lon,lat,degree) #114.4897km
#' 
#' @export degree2Km
`degree2Km` <-
function(lon,lat,degree){
                      x1 <- lon
                      y1 <- lat
                      
                      a <- cos(y1*pi/180)*cos(y1*pi/180)*sin((1*pi/180)/2)*sin((1*pi/180)/2);
                        c <- 2*atan2(sqrt(a),sqrt(1-a));
                        R <- 6371;
                        dx1 <- R*c
                        
              return(dx1 * degree)}

