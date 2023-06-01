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
#' @seealso \code{\link{st_area}}
#' @examples
#' 
#' xrange  <- c(0,4)
#' yrange  <- c(52,56)
#' resx    <- 0.25
#' resy    <- 0.125
#' 
#' #-create grid and assign value column
#' grd <- createGrid(xrange,yrange,resx,resy,type="GridDF",exactBorder=TRUE)
#' st_crs(grd) <- 4326
#' grd$value <- runif(nrow(st_coordinates(st_centroid(grd))),5,10)
#' 
#' #- find gridcells with maximum surface equal or smaller to 'threshold' and return their
#' #  total value and gridcells involved (diagonal means if area may include x+1,y+1 steps v
#' #  versus only x+1 or y+1 steps). Note, the larger the threshold, the longer the function
#' #  will run!
#' res     <- findArea(grd,threshold=1000,diagonal=TRUE)
#' 
#' #- Plot the result
#' plot(st_geometry(grd),type="p",pch=19,cex=0.5)
#' axis(1); axis(2)
#' selec   <- which.min(res$minval)
#' points(st_coordinates(st_centroid(grd))[as.numeric(unlist(strsplit(res[selec,"idxs"]," "))),],col=2,lwd=3)
#' 
#' @export findArea
findArea <- function(grid,threshold=100,diagonal=TRUE){

  if(is.na(st_crs(grid))) stop("CRS is set as NA, threshold is measured in km, so need a CRS to be set")
  if(!"sf" %in% class(grid)) stop(paste("Function not defined for class",class(grid)))
  if(!"value" %in% colnames(grid)) stop("No 'value' column available in data slot")
  grid$surface    <- st_area(grid)/(1000*1000)
  coords          <- st_coordinates(st_centroid(grid))
  
  
  cellsize  <- apply(abs(apply(subset(as.data.frame(st_coordinates(grid)),L2==st_coordinates(grid)[1,"L2"])[,c("X","Y")],2,diff)),2,max,na.rm=T)
  
  resx      <- cellsize["X"]
  resy      <- cellsize["Y"]
  
  storeVals     <- data.frame(minval=rep(0,nrow(coords)),idxs=rep(NA,nrow(coords)))
  for(i in 1:nrow(coords)){
    totSurf     <- 0
    pts         <- i
    while(an(totSurf) <= threshold){
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
      pts       <- unique(c(pts,idx[which.min(grid$value[idx])]))
      totSurf   <- sum(grid$surface[pts])
    }
    storeVals$minval[i] <- grid$value[pts]
    storeVals$idxs[i]   <- paste(pts,collapse= " ")
  }
  return(storeVals)}

