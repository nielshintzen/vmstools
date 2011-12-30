`createGrid` <-
function(xrange
                             ,yrange
                             ,resx
                             ,resy
                             ,type="GridTopology"
                             ){
                
                library(sp)             
                roundDigit  <- max(getndp(resx),getndp(resy),na.rm=T)
                xborder     <- seq(floor(sort(xrange)[1]),ceiling(sort(xrange)[2]),resx)
                yborder     <- seq(floor(sort(yrange)[1]),ceiling(sort(yrange)[2]),resy)
                
                grid        <- GridTopology(c(xborder[1],yborder[1]),c(resx,resy),c(length(xborder),length(yborder)))
                if(type=="SpatialGrid"){
                  grid      <- SpatialGrid(grid=grid);
                }
                if(type=="SpatialPixels"){
                  grid      <- SpatialGrid(grid=grid);
                  gridded(grid) = TRUE
                  grid      <- as(grid,"SpatialPixels");
                }
                if(type=="SpatialPixelsDataFrame"){
                  grid      <- SpatialGrid(grid=grid);
                  gridded(grid) = TRUE
                  grid      <- as(grid,"SpatialPixels");
                  sPDF      <- as(grid,"SpatialPixelsDataFrame")
                  sPDF@data     <- data.frame(rep(0,nrow(coordinates(sPDF))))
                  colnames(sPDF@data) <- "data"
                  grid          <- sPDF
                }
                if(type=="SpatialGridDataFrame"){
                  grid      <- SpatialGrid(grid=grid);
                  sPDF          <- as(grid,"SpatialGridDataFrame")
                  sPDF@data     <- data.frame(rep(0,nrow(coordinates(sPDF))))
                  colnames(sPDF@data) <- "data"
                  grid          <- sPDF
                }
                
              return(grid)}

