`createGrid` <-
function(xrange
                             ,yrange
                             ,resx
                             ,resy
                             ,type="GridTopology"
                             ){
                
                library(sp)             
                if(((xrange[2]-xrange[1])%%resx)==0 && ((yrange[2]-yrange[1])%%resy)==0){
                  roundDigit  <- max(getndp(resx),getndp(resy),na.rm=T)+1  
                }else {
                  roundDigit  <- max(getndp(resx),getndp(resy),na.rm=T)-1
                }
                xborder     <- round(seq(floor(sort(xrange)[1]*(10^roundDigit))/(10^roundDigit),ceiling(sort(xrange)[2]*(10^roundDigit))/(10^roundDigit),resx),roundDigit)
                yborder     <- round(seq(floor(sort(yrange)[1]*(10^roundDigit))/(10^roundDigit),ceiling(sort(yrange)[2]*(10^roundDigit))/(10^roundDigit),resy),roundDigit)
                
                if(xrange[1]<xborder[1] | xrange[2]>rev(xborder)[1] | yrange[1]<yborder[1] | yrange[2]>rev(yborder)[1]) stop("Grid range smaller than specified bounds (bug!)")
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

