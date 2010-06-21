`createGrid` <-
function(xrange
                             ,yrange
                             ,resx
                             ,resy
                             ){
                
                library(sp)             
                roundDigit  <- max(getndp(resx),getndp(resy),na.rm=T)
                xborder     <- seq(floor(sort(xrange)[1]),ceiling(sort(xrange)[2]),resx)
                yborder     <- seq(floor(sort(yrange)[1]),ceiling(sort(yrange)[2]),resy)
                
                grid        <- GridTopology(c(xborder[1],yborder[1]),c(resx,resy),c(length(xborder),length(yborder)))
              return(grid)}

