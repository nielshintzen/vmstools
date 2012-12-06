`createGrid` <-
function(xrange
                             ,yrange
                             ,resx
                             ,resy
                             ,type="GridTopology"
                             ,exactBorder=F
                             ){
                
                require(sp)
                if(exactBorder){
                  roundDigit    <- max(getndp(resx),getndp(resy),na.rm=T)+1
                  xborder       <- round(seq(xrange[1]+resx/2,xrange[2]-resx/2,resx),roundDigit)
                  yborder       <- round(seq(yrange[1]+resy/2,yrange[2]-resy/2,resy),roundDigit)
                } else {
                    roundDigit  <- max(getndp(resx),getndp(resy),na.rm=T)
                    xborder     <- round(seq(floor(sort(xrange)[1]*(10^roundDigit))/(10^roundDigit),ceiling(sort(xrange)[2]/resx)*resx,resx),roundDigit)
                    yborder     <- round(seq(floor(sort(yrange)[1]*(10^roundDigit))/(10^roundDigit),ceiling(sort(yrange)[2]/resy)*resy,resy),roundDigit)
                 }

                if(!exactBorder){
                  if(round(xrange[1],roundDigit)<xborder[1] | round(xrange[2],roundDigit)>rev(xborder)[1] |
                     round(yrange[1],roundDigit)<yborder[1] | round(yrange[2],roundDigit)>rev(yborder)[1]) stop("Grid range smaller than specified bounds (bug!)")
                }
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

