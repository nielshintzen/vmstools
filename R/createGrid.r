`createGrid` <-
function(xrange
                             ,yrange
                             ,resx
                             ,resy
                             ,type="GridTopology"
                             ,exactBorder=FALSE
                             ){
                
                require(sp)
                if(exactBorder){
                  xs            <- seq(xrange[1]+resx/2,xrange[2]+resx,resx)
                  xborder       <- range(xs)
                  ys            <- seq(yrange[1]+resy/2,yrange[2]+resy,resy)
                  yborder       <- range(ys)
                  } else {
                    roundDigitx <- getndp(resx)
                    roundDigity <- getndp(resy)
                    xs          <- seq(xrange[1],xrange[2] +resx,resx)
                    xborder     <- c((min(floor(xs*10^roundDigitx))/10^roundDigitx)- resx/2,(max(ceiling(xs*10^roundDigitx)) / 10^roundDigitx)+resx/2)
                    ys          <- seq(yrange[1],yrange[2] +resy,resy)
                    yborder     <- c((min(floor(ys*10^roundDigity))/10^roundDigity)- resy/2,(max(ceiling(ys*10^roundDigity)) / 10^roundDigity)+resy/2)
                   }

                grid        <- GridTopology(c(xborder[1],yborder[1]),c(resx,resy),c(length(seq(xborder[1],xborder[2],resx)),length(seq(yborder[1],yborder[2],resy))))
                if(type=="SpatialGrid"){
                  grid      <- SpatialGrid(grid=grid);
                }
                if(type=="SpatialPixels"){
                  grid      <- GridTopology(c(xborder[1]-resx/2,yborder[1]-resy/2),c(resx,resy),c(length(seq(xborder[1]-resx/2,xborder[2]+resx/2,resx)),length(seq(yborder[1]-resy/2,yborder[2]+resy/2,resy))))
                  grid      <- SpatialGrid(grid=grid);
                  gridded(grid) = TRUE
                  grid      <- as(grid,"SpatialPixels");
                }
                if(type=="SpatialPixelsDataFrame"){
                  grid      <- GridTopology(c(xborder[1]-resx/2,yborder[1]-resy/2),c(resx,resy),c(length(seq(xborder[1]-resx/2,xborder[2]+resx/2,resx)),length(seq(yborder[1]-resy/2,yborder[2]+resy/2,resy))))
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
                
               if(!class(grid)=="GridTopology"){
                 bb <- bbox(grid)
                 if(bb[1,1] > min(xrange) | max(xrange) > bb[1,2] | bb[2,1] > min(yrange) | max(yrange) > bb[2,2])
                  stop("Dimensions of grid not large enough to span xranges and yranges (bug)")
               }
                
              return(grid)}

