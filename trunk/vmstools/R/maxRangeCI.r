

maxRangeCI <- function(x,y,Time,speed){

  #Pre-Calculation to speed up the code
  pi180 <- pi/180
  cosy1 <- cos(y[1]*pi180)
  cosy2 <- cos(y[2]*pi180)
  
  #Calculate maximum distance in km
  dmax <- Time/60*sum(speed,na.rm=TRUE)/2*1.852
  
  #Calculate d from Haversine function
  d   <- distance(x[1],y[1],x[2],y[2])
  
  #Calculate a and b as in Mills et al. 2006 paper
    warn<- 0
  if(d >= dmax){
    warning(paste("Distance too far to reach with current speed estimate ",round(x,3)," ",round(y,3),"\n"))
    dmax <- d
    warn <- 1
  }
  a <- dmax/2
  b <- sqrt((dmax^2 - d^2)/4)
  
  if(d == 0){
    o <- 0
  } else {
    dx      <- (x[2] - x[1])*pi180
    dy      <- (y[2] - y[1])*pi180
    o       <- atan2(sin(dx)*cosy2,cosy1*sin(y[2]*pi180)-sin(y[1]*pi180)*cosy2*cos(dx))
    angles  <- (o*(180/pi)) %% 360
    
    angle2  <- ifelse(angles >= 0 & angles < 180, 90 - angles,270-angles)
    o       <- angle2*(pi180)
  }
  #See also: http://www.movable-type.co.uk/scripts/latlong.html
  Bx        <- cosy2*cos((x[2]-x[1])*pi180)
  By        <- cosy2*sin((x[2]-x[1])*pi180)
  mid.x     <- (x[1]*pi180) + atan2(By,cosy1+Bx)
  mid.y     <- atan2(sin(y[1]*pi180) + sin(y[2]*pi180),sqrt((cosy1+Bx)^2 + By^2))
  mid.x     <- mid.x*180/pi
  mid.y     <- mid.y*180/pi

  a         <- c(km2Degree(mid.x,mid.y,a),a/111.2)
  b         <- c(km2Degree(mid.x,mid.y,b),b/111.2)
  
  #See also Pfoser and Jensen 1999 Capturing the Uncertainty of Moving-Object representation
  u         <- 0:360*pi180
  xres      <- mid.x + a[1] * cos(o) * cos(u) - b[1] * sin(o) * sin(u)
  yres      <- mid.y + a[2] * sin(o) * cos(u) + b[2] * cos(o) * sin(u)

  return(list(matrix(c(xres,yres),ncol=2),dmax,warn))}