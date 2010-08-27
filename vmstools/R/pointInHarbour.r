pointInHarbour <- function(lon,lat,harbours){

    inHarbour   <- numeric()
    for(hars in 1:length(harbours)){
      print(hars)
      lonDegree <- 0.1*(harbours[[hars]]$range/(11.2*lonLatRatio(lon=harbours[[hars]]$lon,lat=harbours[[hars]]$lat)[1]))
      latDegree <- 0.1*(harbours[[hars]]$range/11.2)

      #WHY IS LEVEL ~0.392??
      harPolygon <- ellipse(x=0,scale=c(lonDegree,latDegree),centre=c(harbours[[hars]]$lon,harbours[[hars]]$lat),level=0.392)
      res        <- point.in.polygon(point.x=lon,point.y=lat,pol.x=harPolygon[,"x"],pol.y=harPolygon[,"y"])
      inHarbour  <- unique(c(inHarbour,which(res == 1)))
    }

return(inHarbour)}
    
