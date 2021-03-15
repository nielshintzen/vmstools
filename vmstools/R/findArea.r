#' Find area surface smaller than threshold
#' 
#' Find gridcells with maximum surface equal or smaller to 'threshold' and
#' return their total value and gridcells involved
#' 
#' Diagonal means if area may include x+1,y+1 steps v versus only x+1 or y+1
#' steps
#' 
#' @param grid A SpatialGridDataFrame
#' @param threshold Maximum value of surfaces added up
#' @param diagonal Allow diagonal steps in gridcell selection
#' @return data.frame with minimum surface area and gridcells selected
#' @author Niels T. Hintzen
#' @seealso \code{\link{surface}}
#' @examples
#' 
#' xrange  <- c(0,4)
#' yrange  <- c(52,56)
#' resx    <- 0.25
#' resy    <- 0.125
#' 
#' #-create grid and assign value column
#' grd <- createGrid(xrange,yrange,resx,resy,type="SpatialGridDataFrame",exactBorder=TRUE)
#' grd@data$value <- runif(nrow(coordinates(grd)),5,10)
#' 
#' #- find gridcells with maximum surface equal or smaller to 'threshold' and return their
#' #  total value and gridcells involved (diagonal means if area may include x+1,y+1 steps v
#' #  versus only x+1 or y+1 steps). Note, the larger the threshold, the longer the function
#' #  will run!
#' res     <- findArea(grd,threshold=1000,diagonal=TRUE)
#' 
#' #- Plot the result
#' plot(grd,type="p",pch=19,cex=0.5)
#' map.axes()
#' selec   <- which.min(res$minval)
#' points(coordinates(grd)[as.numeric(unlist(strsplit(res[selec,"idxs"]," "))),],col=2,lwd=3)
#' 
#' @export findArea
findArea <- function(grid,threshold=100,diagonal=TRUE){

  if(class(grid) != "SpatialGridDataFrame") stop(paste("Function not defined for class",class(grid)))
  if(!"value" %in% colnames(grid@data)) stop("No 'value' column available in data slot")
  grid    <- surface(grid)
  coords  <- coordinates(grid)
  
  resx    <- grid@grid@cellsize[1]
  resy    <- grid@grid@cellsize[2]
  
  storeVals     <- data.frame(minval=rep(0,nrow(coords)),idxs=rep(NA,nrow(coords)))
  for(i in 1:nrow(coords)){
    totSurf     <- 0
    pts         <- i
    while(totSurf <= threshold){
      dists     <- do.call(pmin, lapply(as.list(pts), function(x) {
                    distance(coords[x, 1], coords[x, 2], coords[-pts, 1], coords[-pts,2])}))
      if(diagonal==TRUE) idx       <- (1:nrow(coords))[-pts][which(dists <= do.call(max, lapply(as.list(pts), function(x) {
                                                max(distance(coords[x, 1], coords[x, 2], coords[x, 1]-resx, coords[x,2]-resy),
                                                    distance(coords[x, 1], coords[x, 2], coords[x, 1]+resx, coords[x,2]+resy))})))]
      if(diagonal==FALSE) idx       <- (1:nrow(coords))[-pts][which(dists <= do.call(max, lapply(as.list(pts), function(x) {
                                                max(distance(coords[x, 1], coords[x, 2], coords[x, 1]-resx, coords[x,2]),
                                                    distance(coords[x, 1], coords[x, 2], coords[x, 1]+resx, coords[x,2]),
                                                    distance(coords[x, 1], coords[x, 2], coords[x, 1], coords[x,2]-resy),
                                                    distance(coords[x, 1], coords[x, 2], coords[x, 1], coords[x,2]+resy))})))]
      pts       <- unique(c(pts,idx[which.min(grid@data$value[idx])]))
      totSurf   <- sum(grid@data$cellArea[pts])
    }
    storeVals$minval[i] <- grid@data$value[pts]
    storeVals$idxs[i]   <- paste(pts,collapse= " ")
  }
  return(storeVals)}
