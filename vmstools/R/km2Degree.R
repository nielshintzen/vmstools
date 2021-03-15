#' Compute distance from kilometers into degrees
#' 
#' Function transformes the distance expressed in kilometers into degrees. This
#' based on the GPS location of a point.
#' 
#' 
#' @param lon Longitude of the GPS position
#' @param lat Latitude of the GPS positiona
#' @param km Value in Km to turn into degrees
#' @note Computation of degrees is approximation based on the Haversine formula
#' @author Niels T. Hintzen
#' @seealso \code{\link{distance}}, \code{\link{degree2Km}},
#' \code{\link{lonLatRatio}}
#' @references EU lot 2 project
#' @examples
#' 
#' lon <- -4
#' lat <- 50
#' km  <- 114.4897
#' 
#' km2Degree(lon,lat,km) #1.601833
#' 
#' @export km2Degree
`km2Degree` <-
function(lon,lat,km){
                      x1 <- lon
                      y1 <- lat
                      
                      a <- cos(y1*pi/180)*cos(y1*pi/180)*sin((1*pi/180)/2)*sin((1*pi/180)/2);
                        c <- 2*atan2(sqrt(a),sqrt(1-a));
                        R <- 6371;
                        dx1 <- R*c
                        
              return(km / dx1)}

