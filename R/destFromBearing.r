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