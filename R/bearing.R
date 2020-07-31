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