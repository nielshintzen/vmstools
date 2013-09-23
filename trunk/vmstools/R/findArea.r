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
