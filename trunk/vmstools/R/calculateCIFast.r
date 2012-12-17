calculateCIFast <- function(    int
                               ,tacint
                               ,params
                               ,grid
                               ,plot=F){
                               
  mxr     <-  maxRangeCI(x  =c(int[-1,1][1],rev(int[-1,1])[1]),
                         y  =c(int[-1,2][1],rev(int[-1,2])[1]),
                         time.=c(difftime(tacint$SI_DATIM[2],tacint$SI_DATIM[1],units="mins")),
                         speed=pmax(tacint$SI_SP,rep(distanceInterpolation(list(int)) / 1.852 /
                                                     c(difftime(tacint$SI_DATIM[2],tacint$SI_DATIM[1],units="hours")),2)))
                                                     
  if(plot){
    par(mfrow=c(2,2))
    plot(mxr[[1]][,1],mxr[[1]][,2],type="l",xlab="Longitude",ylab="Latitude",asp=1/lonLatRatio(mxr[[1]][1,1],mxr[[1]][1,2]),main=paste(int[1,]))
    lines(int[-1,1],int[-1,2],col=2)
  }
  
  xrange <- range(mxr[[1]][,1]); yrange <- range(mxr[[1]][,2])
  xrg    <- range(int[-1,1]); yrg <- range(int[-1,2])
  if(xrange[1] > xrg[1]) xrange[1] <- xrg[1] - diff(xrg)*0.1
  if(xrange[2] < xrg[2]) xrange[2] <- xrg[2] + diff(xrg)*0.1
  if(yrange[1] > yrg[1]) yrange[1] <- yrg[1] - diff(yrg)*0.1
  if(yrange[2] < yrg[2]) yrange[2] <- yrg[2] - diff(yrg)*0.1

  grid      <- createGrid(c((xrange[1] - grid@cellcentre.offset[1])                   %/%grid@cellsize[1] * grid@cellsize[1] + grid@cellcentre.offset[1],
                            (xrange[2] - grid@cellcentre.offset[1] + grid@cellsize[1])%/%grid@cellsize[1] * grid@cellsize[1] + grid@cellcentre.offset[1]),
                          c((yrange[1] - grid@cellcentre.offset[2])                   %/%grid@cellsize[2] * grid@cellsize[2] + grid@cellcentre.offset[2],
                            (yrange[2] - grid@cellcentre.offset[2] + grid@cellsize[2])%/%grid@cellsize[2] * grid@cellsize[2] + grid@cellcentre.offset[2]),
                          grid@cellsize[1],grid@cellsize[2],type="SpatialGridDataFrame")
  coords    <- coordinates(grid)

  # Distance to begin or end point
  bpDistan  <- distance(tacint$SI_LONG[1],tacint$SI_LATI[1],coords[,1],coords[,2])
  epDistan  <- distance(tacint$SI_LONG[2],tacint$SI_LATI[2],coords[,1],coords[,2])
  pdistan   <- pmin(bpDistan,epDistan)
  
  if(plot){ image(t(matrix(pdistan,ncol=grid@grid@cells.dim[1],nrow=grid@grid@cells.dim[2],byrow=T)[grid@grid@cells.dim[2]:1,]),col=rev(heat.colors(12))); box()}

  # Distance to interpolation
  intDistan <- do.call(pmin,lapply(as.list(2:nrow(int)),function(x){distance(int[x,1],int[x,2],coords[,1],coords[,2])}))
    
  if(plot){ image(t(matrix(intDistan,ncol=grid@grid@cells.dim[1],nrow=grid@grid@cells.dim[2],byrow=T)[grid@grid@cells.dim[2]:1,]),col=rev(heat.colors(12))); box()}

  CI        <- N1p0(intDistan*params$distscale,0,pdistan^params$sigline,0)
  if(max(CI,na.rm=T) < 0.1) warning("Prediction max(CI) is very small")

  if(plot){ image(t(matrix(CI,ncol=grid@grid@cells.dim[1],nrow=grid@grid@cells.dim[2],byrow=T)[grid@grid@cells.dim[2]:1,]),col=rev(heat.colors(12))); box()}

  grid@data$data <- CI
return(grid)}
                               
