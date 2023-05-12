#' Calculate the CI between two succeeding VMS datapoints
#' 
#' The interpolated tracks can be surrounded by a sort of confidence interval
#' representing the outer region a vessel could have travelled between two
#' succeeding datapoints. Within this function the CI's are computed.
#' 
#' 
#' @param int interpolation, as data.frame from 'interpolateTacsat'
#' @param tacint tacsat records (two rows) corresponding with the interpolation
#' @param params list of parameters used to perform interpolation
#' @param grid object of class 'GridTopology' specifying the grid dimensions
#' @param plot Logical. Whether the result of the interpolation must be plotted
#' @return Returns the Confidence Interval on a grid of class
#' 'SpatialGridDataFrame' with the CI values in the data slot.
#' @note Computation can take a very long time if either grid resolution is
#' high or if many interpolations are used.
#' @author Niels T. Hintzen
#' @seealso \code{\link{interpolateTacsat}},\code{\link{maxRangeCI}}
#' @references Hintzen et al. 2010 Improved estimation of trawling tracks using
#' cubic Hermite spline interpolation of position registration data, EU lot 2
#' project
#' @examples
#' 
#' data(tacsat)
#' 
#'   #Sort the Tacsat data
#' tacsat     <- sortTacsat(tacsat)
#' tacsat     <- tacsat[1:1000,]
#' 
#'   #Filter the Tacsat data
#' tacsat          <- filterTacsat(tacsat,c(2,6),hd=NULL,remDup=TRUE)
#' 
#'   #Interpolate the VMS data
#' interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
#'                     res=100,method="cHs",params=list(fm=0.5,distscale=20,
#'                     sigline=0.2,st=c(2,6)),headingAdjustment=0)
#' 
#'   #Create the final grid where all interpolations should fit on
#' xrange        <- c(2,3); yrange <- c(51,52)
#' grid          <- createGrid(xrange,yrange,resx=0.01,resy=0.005)
#' 
#' res           <- calculateCI(interpolation[[4]],
#'                              tacsat[interpolation[[4]][1,],],
#'                              params=list(fm=0.25,distscale=3.1,sigline=0.4,st=c(2,6)),
#'                              grid=grid,
#'                              plot=TRUE)
#' 
#' @export calculateCI
calculateCI <- function(    int
                               ,tacint
                               ,params
                               ,grid
                               ,spatialGrid
                               ,plot=FALSE){

  if (!"SI_DATIM" %in% colnames(tacint))
        tacint$SI_DATIMIM <- as.POSIXct(paste(tacint$SI_DATE, tacint$SI_TIME,
            sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M")
  mxr     <-  maxRangeCI(x  =c(int[-1,1][1],rev(int[-1,1])[1]),
                         y  =c(int[-1,2][1],rev(int[-1,2])[1]),
                         Time=an(difftime(tacint$SI_DATIM[2],tacint$SI_DATIM[1],units="mins")),
                         speed=pmax(tacint$SI_SP,rep(distanceInterpolation(list(int)) / 1.852 /
                                                     an(difftime(tacint$SI_DATIM[2],tacint$SI_DATIM[1],units="hours")),2)))
                                                     
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

  
  cellsize  <- apply(abs(apply(subset(as.data.frame(st_coordinates(grid)),L2==st_coordinates(grid)[1,"L2"])[,c("X","Y")],2,diff)),2,max,na.rm=T)
  celldim   <- c(length(apply(st_coordinates(grid),2,unique)$X),length(apply(st_coordinates(grid),2,unique)$Y))
  celloffset<- st_bbox(grid)[c("xmin","ymin")] - cellsize/2

  newxrange <- c((xrange[1] - celloffset[1])                   %/%cellsize[1] * cellsize[1] + celloffset[1],
                            (xrange[2] - celloffset[1] + cellsize[1])%/%cellsize[1] * cellsize[1] + celloffset[1])
  newyrange <- c((yrange[1] - celloffset[2])                   %/%cellsize[2] * cellsize[2] + celloffset[2],
                            (yrange[2] - celloffset[2] + cellsize[2])%/%cellsize[2] * cellsize[2] + celloffset[2])

  origxs    <- seq(celloffset[1],celloffset[1] + cellsize[1] * (celldim[1]),by= cellsize[1])
  origys    <- seq(celloffset[2],celloffset[2] + cellsize[2] * (celldim[2]),by= cellsize[2])

  if((min(newxrange) > st_bbox(grid)["xmax"]) | (max(newxrange) < st_bbox(grid)["xmin"])) stop("CI outside grid")
  if((min(newyrange) > st_bbox(grid)["ymax"]) | (max(newyrange) < st_bbox(grid)["ymin"])) stop("CI outside grid")

  subxs     <- apply(abs(outer(origxs,newxrange,"-")),2,which.min)
  subys     <- apply(abs(outer(origys,newyrange,"-")),2,which.min)

  idxy      <- expand.grid(x=1:(celldim[1]-1),y=1:(celldim[2]-1)) 
  gridsm    <- grid[which(idxy[,"x"] %in% seq(subxs[1],subxs[2]) & idxy[,"y"] %in% seq(subys[1],subys[2]))]
  coords    <- st_coordinates(st_centroid(gridsm))

  # Distance to begin or end point
  bpDistan  <- distance(tacint$SI_LONG[1],tacint$SI_LATI[1],coords[,1],coords[,2])
  epDistan  <- distance(tacint$SI_LONG[2],tacint$SI_LATI[2],coords[,1],coords[,2])
  pdistan   <- pmin(bpDistan,epDistan)
  
  if(plot){ image(t(matrix(pdistan,ncol=length(seq(subxs[1],subxs[2])),nrow=length(seq(subys[1],subys[2])),byrow=TRUE)),col=rev(heat.colors(12))); box()}

  # Distance to interpolation
  intDistan <- do.call(pmin,lapply(as.list(2:nrow(int)),function(x){distance(int[x,1],int[x,2],coords[,1],coords[,2])}))
    
  if(plot){ image(t(matrix(intDistan,ncol=length(seq(subxs[1],subxs[2])),nrow=length(seq(subys[1],subys[2])),byrow=TRUE)),col=rev(heat.colors(12))); box()}

  CI        <- N1p0(intDistan*params$distscale,0,pdistan^params$sigline,0)
  if(max(CI,na.rm=TRUE) < 0.1) warning("Prediction max(CI) is very small")

  if(plot){ image(t(matrix(CI,ncol=length(seq(subxs[1],subxs[2])),nrow=length(seq(subys[1],subys[2])),byrow=TRUE)),col=rev(heat.colors(12))); box()}
  gridsm <- st_as_sf(data=CI,gridsm)  
return(gridsm)}

