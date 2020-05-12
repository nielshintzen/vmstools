ICESrectangle2CSquare <- function(rectangles,degrees,onLand=T){

  lonlats     <- ICESrectangle2LonLat(rectangles)
  
  if(degrees >= 1)
    ret       <- CSquare(lonlats$SI_LONG,lonlats$SI_LATI,degrees=degrees)
  if(degrees < 1){
    squares   <- lapply(as.list(1:nrow(lonlats)),function(x){
                                                 return(expand.grid(SI_LONG=seq(lonlats[x,"SI_LONG"],lonlats[x,"SI_LONG"]+1-degrees,degrees)+degrees/2,SI_LATI=seq(lonlats[x,"SI_LATI"],lonlats[x,"SI_LATI"]+0.5-degrees,degrees)+degrees/2))})
    ret       <- lapply(squares,function(x){return(CSquare(x$SI_LONG,x$SI_LATI,degrees=degrees))})
  }
  names(ret)  <- rectangles
  if(onLand==F){
    if(!exists("ICESareas"))
      data(ICESareas)
    require(sp)
    idx       <- lapply(ret,function(x){idx <- over(SpatialPoints(CSquare2LonLat(x,degrees)[,c("SI_LONG","SI_LATI")]-1e-05),ICESareas);
                            return(which(is.na(idx)==F))})  #1e-5 is rounding error introduced in CSquare2LonLat
    res       <- lapply(as.list(1:length(ret)),function(x){return(ret[[x]][idx[[x]]])})
    names(res)<- names(which(unlist(lapply(idx,length))>0))
    ret       <- res
  }
  return(ret)}

#- Example
#degrees       <- 0.05
#rectangles    <- c("33F3","33F4","32F3","34F3")
#squares       <- ICESrectangle2CSquare(rectangles,degrees)
#squares       <- do.call(c,squares)
#par(xaxs="i",yaxs="i",las=1,oma=c(2,2,1,1))
#plot(ICESareas,xlim=c(2,5),ylim=c(51.5,53),col="lightblue",fill=T);
#map.axes()
#abline(v=seq(2,5,1),lty=3)
#abline(h=seq(51.5,53,0.5),lty=3)
#lonLatSquares <- CSquare2LonLat(squares,degrees=degrees)
#points(lonLatSquares$SI_LONG,lonLatSquares$SI_LATI,pch=19,cex=0.2,col=2)
#
#squares       <- ICESrectangle2CSquare(rectangles,degrees,onLand=F)
#squares       <- do.call(c,squares)
#lonLatSquares <- CSquare2LonLat(squares,degrees=degrees)
#points(lonLatSquares$SI_LONG,lonLatSquares$SI_LATI,pch=19,cex=0.2,col=4)

