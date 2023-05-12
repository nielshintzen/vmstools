#' Compute bearing between two points on a sphere (approximation of the earth)
#' at the starting point
#' 
#' Compute the bearing between two GPS locations defined in longitude and
#' latitude notation on the earth. The earth is assumed to have a perfect
#' spherical shape. Bearing is returned in compass degrees.
#' 
#' 
#' @param lon Longitude of point 2
#' @param lat Latitude of point 2
#' @param lonRef Longitude of point 1
#' @param latRef Latitude of point 1
#' @author Niels T. Hintzen
#' @seealso \code{\link{km2Degree}}, \code{\link{degree2Km}},
#' \code{\link{lonLatRatio}},\code{\link{distance}}
#' @references EU Lot 2 project, based on the Haversine formula, see also:
#' Hintzen et al. 2010 Fisheries Research
#' @examples
#' 
#' lon <- -4
#' lat <- 50
#' lonRef <- -4.2
#' latRef <- 51
#' 
#' bearing(lon,lat,lonRef,latRef) #352.8271
#' 
#' @export bearing
bearing <- function(lon,lat,lonRef,latRef){

                    x1  <- lon
                    y1  <- lat
                    x2  <- lonRef
                    y2  <- latRef
                    
                    y   <- sin((x2-x1)*pi/180) * cos(y2*pi/180)
                    x   <- cos(y1*pi/180) * sin(y2*pi/180) - sin(y1*pi/180) * cos(y2*pi/180) * cos((x2-x1)*pi/180)
                    bearing <- atan2(y,x)*180/pi
                    bearing <- (bearing + 360)%%360
          return(bearing)}

                    
#bearing(3.15,51.64,3.15,51.65)

#hello world


lat <- c(53.63125, 53.62593, 53.62003, 53.61626, 53.61683, 53.61750)
lon <- c(-6.17412, -6.17330, -6.17629, -6.17105, -6.17068, -6.17085)
plot(lon,lat,pch=19,col=2)
points(destFromBearing(lon[-length(lon)],lat[-length(lat)],bearing(lon[-length(lon)],lat[-length(lat)],lon[-1],lat[-1]),distance=0.1),pch=1)
text(x = lon, y = lat,
     labels = as.character(1:length(lat)), adj = -0.8)


bearing(lon[1],lat[1],lon[2],lat[2])
