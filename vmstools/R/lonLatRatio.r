#' Calculate the ratio between 1 degree in longitude versus 1 degree in
#' latitude
#' 
#' The distance in Km on the longitude direction changes along the latitude
#' direction. This function computes the ratio between 1 degree in the
#' longitude direction depending on the latitude of the GPS position. Returns
#' the ratio's of two GPS locations (two succeeding VMS datapoints). Can be
#' used with 1 GPS position too, return NA for second value.
#' 
#' 
#' @param lon Longitude of the two GPS positions
#' @param lat Latitude of the two GPS positions
#' @note Computation is approximation based on the Haversine formula
#' @author Niels T. Hintzen
#' @seealso \code{\link{distance}}, \code{\link{degree2Km}},
#' \code{\link{km2Degree}}
#' @references EU lot 2 project
#' @examples
#' 
#' lon <- -4
#' lat <- 50
#' 
#' lonLatRatio(lon,lat)
#' 
#' @export lonLatRatio
`lonLatRatio` <-
    function(x1,lat){
      #Based on the Haversine formula
      #At the position, the y-position remains the same, hence, cos(lat)*cos(lat) instead of cos(lat) * cos(y2)
      a <- cos(lat*pi/180)*cos(lat*pi/180)*sin((0.1*pi/180)/2)*sin((0.1*pi/180)/2);
      c <- 2*atan2(sqrt(a),sqrt(1-a));
      R <- 6371;
      dx1 <- R*c

    return(c(dx1/11.12))}
