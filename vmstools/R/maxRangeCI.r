maxRangeCI <- function(lon,lat,time.,speed){

                    x1 <- lon[1]
                    x2 <- lon[2]
                    y1 <- lat[1]
                    y2 <- lat[2]

                      #Calculate maximum distance in km
                    dmax    <- time./60*mean(speed,na.rm=T)*1.852

                      #Calculate d from Haversine function
                    aH  <- sin(((y2-y1)*pi/180)/2)*sin(((y2-y1)*pi/180)/2) + cos(y1*pi/180) * cos(y2*pi/180) *
                           sin(((x2-x1)*pi/180)/2) * sin(((x2-x1)*pi/180)/2)
                    c   <- 2*atan2(sqrt(aH),sqrt(1-aH))
                    d   <- 6371*c
                    
                      #Calculate a and b as in Mills et al. 2006 paper
                    a   <- dmax/2
                    warn<- 0
                    if(d >= dmax){
                      warning(paste("Distance too far to reach with current speed estimate ",round(lon,3)," ",round(lat,3),"\n"))
                      dmax <- d
                      warn <- 1
                    }
                    b   <- sqrt((dmax^2 - d^2)/4)
                    
                    if(d == 0){
                      o <- 0
                    } else {
                        dx      <- (x2 - x1)*pi/180
                        dy      <- (y2 - y1)*pi/180
                        o       <- atan2(sin(dx)*cos(y2*pi/180),cos(y1*pi/180)*sin(y2*pi/180)-sin(y1*pi/180)*cos(y2*pi/180)*cos(dx))
                        angles  <- o*(180/pi)

                        angles <- (angles + 360)%%360

                        if(angles >=0   & angles < 90)  angle2 <- 90-angles
                        if(angles >=90  & angles < 180) angle2 <- angles - 90
                        if(angles >=180 & angles < 270) angle2 <- 90 - (angles-180)
                        if(angles >=270 & angles < 360) angle2 <- (angles - 270)

                        o <- angle2*(pi/180)
                    }
                      #See also: http://www.movable-type.co.uk/scripts/latlong.html
                    Bx    <- cos(y2*pi/180)*cos((x2-x1)*pi/180)
                    By    <- cos(y2*pi/180)*sin((x2-x1)*pi/180)
                    mid.x <- (x1*pi/180) + atan2(By,cos(y1*pi/180)+Bx)
                    mid.y <- atan2(sin(y1*pi/180) + sin(y2*pi/180),sqrt((cos(y1*pi/180)+Bx)^2 + By^2))
                    mid.x <- mid.x*180/pi
                    mid.y <- mid.y*180/pi
                    x <- numeric()
                    y <- numeric()

                    a <- c(km2Degree(mid.x,mid.y,a),a/111.2)
                    b <- c(km2Degree(mid.x,mid.y,b),b/111.2)

                      #See also Pfoser and Jensen 1999 Capturing the Uncertainty of Moving-Object representation
                    for (k in 1:360){
                      u <- k*pi/180
                      x[k] <- mid.x + a[1] * cos(o) * cos(u) - b[1] * sin(o) * sin(u)
                      y[k] <- mid.y + a[2] * sin(o) * cos(u) + b[2] * cos(o) * sin(u)
                    }
                return(list(matrix(c(x,y),ncol=2),dmax,warn))}

