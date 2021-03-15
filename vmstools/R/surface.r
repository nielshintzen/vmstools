#' Calculate surface from grid cells or polygons
#' 
#' Calculate the surface in km2 of the grid / polygons that has been used
#' 
#' Method UTM might take longer due to way of calculation, but is more precise
#' than the Trapezoid function, especially when larger gridcells are used.
#' 
#' @param obj defined SpatialGridDataFrame or SpatialPolygons object (see 'sp'
#' package)
#' @param method Method to be used to calculate surface, either Trapezoid or
#' UTM
#' @param includeNA Whether to include cells which do not hold any data
#' @param zone Include UTM zone notation (or detected automatically when NULL
#' @return If obj is a SpatialGridDataFrame an additional column named
#' 'cellArea' is returned which holds the km2 area of the grid cell. If obj is
#' a SpatialPolygons, each polygon area slot is filled with the area in km2.
#' @author Niels T. Hintzen
#' @seealso \code{\link{createGrid}}
#' @examples
#' 
#' data(tacsat)
#' 
#'   #Sort the Tacsat data
#' tacsat     <- sortTacsat(tacsat)
#' tacsat     <- tacsat[1:1000,]
#' 
#'   #Get the ranges of the tacsat data
#' xrange  <- range(tacsat$SI_LONG,na.rm=TRUE);
#' xrange  <- c(min(xrange) - min(xrange)*0.05,
#'              max(xrange) + max(xrange)*0.05)
#' yrange  <- range(tacsat$SI_LATI,na.rm=TRUE);
#' yrange  <- c(min(yrange) - min(yrange)*0.05,
#'              max(yrange) + max(yrange)*0.05)
#'   #Setup a grid
#' sPDF    <- createGrid(xrange,yrange,resx=0.1,resy=0.05,type="SpatialGridDataFrame")
#' 
#'   #Setup a polygon
#' sP      <- lonLat2SpatialPolygons(lst=list(data.frame(
#'               SI_LONG=c(4,4.5,4.7,4),
#'               SI_LATI=c(54,54,55.5,55.7))))
#' 
#'   #Calculate the cell surface
#' result        <- surface(sPDF,method="Trapezoid",includeNA=TRUE)
#' print(head(result@data))
#' result        <- surface(sP,zone=31)
#' print(result@polygons[[1]]@Polygons[[1]]@area)
#' 
#' @export surface
surface <- function(obj,method="Trapezoid",includeNA=TRUE,zone=NULL){  #Methods can be "Trapezoid" and "UTM"
        require(sp)
        require(PBSmapping)
        if(!class(obj) %in% c('SpatialGridDataFrame','SpatialPolygons'))
          stop("class of obj should be SpatialGridDataFrame or SpatialPolygons")
        
        if(class(obj) == 'SpatialPolygons')
          {
            allSourcePoly <- numeric()
            counter       <- 0
            for(iPol1 in 1:length(obj@polygons)){
              for(iPol2 in 1:length(obj@polygons[[iPol1]])){
                counter   <- counter + 1
                sourcePoly <- data.frame(cbind(1,1:nrow(obj@polygons[[iPol1]]@Polygons[[iPol2]]@coords),
                                         obj@polygons[[iPol1]]@Polygons[[iPol2]]@coords[,1],obj@polygons[[iPol1]]@Polygons[[iPol2]]@coords[,2]))
                rownames(sourcePoly)<-1:nrow(sourcePoly)
                colnames(sourcePoly)<-c("PID","POS","X","Y")
                sourcePoly$PID[]    <-counter
                
                allSourcePoly       <- rbind(allSourcePoly,sourcePoly)
              }
            }
            areas   <- calcArea(as.PolySet(allSourcePoly, projection="LL",zone=zone))$area
            areas   <- data.frame(
                        cbind(areas,do.call(rbind,lapply(obj@polygons,function(x){return(x@labpt)})),
                                    do.call(rbind,lapply(obj@polygons,function(x){return(x@ID)}))),stringsAsFactors=FALSE)
            colnames(areas) <- c("areas","labptx","labpty","ID")


            obj     <- SpatialPolygons(lapply(obj@polygons,function(x){
                                res <- lapply( x@Polygons,function(y){
                                                            subAreas <- subset(areas,ID == x@ID & labptx == ac(y@labpt[1]) & labpty == ac(y@labpt[2]))
                                                            y@area   <- anf(subAreas$areas);return(y)});
                                                          return(Polygons(res,ID=x@ID))}))

          }
        if(!method %in% c("Trapezoid","UTM")) stop("method not available")
        if (class(obj) %in% c('SpatialGridDataFrame')) # not empty...
          {
             if(method == "Trapezoid"){
               res <- ceiling(max(obj@grid@cellsize,na.rm=TRUE)/0.1 * 10)  #automatic scaling
               if(res < 3) res <- 3
               griddims <- summary(obj)$grid
               bboxdims <- bbox(obj)
               stlon    <- bboxdims[1,1]
               stlat    <- bboxdims[2,1]
               enlon    <- bboxdims[1,2]
               enlat    <- bboxdims[2,2]
               sizelon  <- griddims[1,2]
               sizelat  <- griddims[2,2]

               lons     <- seq(stlon,enlon,sizelon)
               lats     <- seq(stlat,enlat,sizelat)

               heights  <- distance(lon=0,lat=stlat,lonRef=0,latRef=stlat+sizelat/(res-1))
               seqlats  <- mapply(seq,lats[1:(length(lats)-1)],lats[2:length(lats)],length.out=res)

               base     <- matrix(mapply(distance,lon=0,lat=c(seqlats),lonRef=sizelon,latRef=c(seqlats)),ncol=res,byrow=TRUE)
               if(dim(base)[1] == 1){
                  base1     <- base[1:(res-1)]
                  base2     <- base[2:res]
                  surface   <- rep(sum(heights * (base1 + base2)/2),each=length(seq(stlon,enlon-sizelon,sizelon)))
               } else {
                  base1    <- base[,1:(res-1)]
                  base2    <- base[,2:res]
                  surface <- rep(apply(heights * (base1 + base2) / 2,1,sum),each=length(seq(stlon,enlon-sizelon,sizelon)))
                }

               obj@data$cellArea <- rev(surface)
            }
          if(method == "UTM"){
            require(PBSmapping)
            griddims <- summary(obj)$grid
            sizelon  <- griddims[1,2]
            sizelat  <- griddims[2,2]

            ltCentreCell<-coordinates(obj)

            for (x in 1:(length(ltCentreCell)/2)){
              if (includeNA) {        # speed up the calculation by dropping cells with fishing=NA  /!| only work for DCF5 and DCF6!
                minX<-ltCentreCell[x,1]-sizelon/2
                maxX<-ltCentreCell[x,1]+sizelon/2
                minY<-ltCentreCell[x,2]-sizelat/2
                maxY<-ltCentreCell[x,2]+sizelat/2
                ltX<-c(minX,minX,maxX,maxX)
                ltY<-c(minY,maxY,maxY,minY)
                sourcePoly<-cbind(rep(1,4),seq(1,4),ltX,ltY)
                rownames(sourcePoly)<-seq(1,4)
                colnames(sourcePoly)<-c("PID","POS","X","Y")

                polyArea<-calcArea(as.PolySet(sourcePoly, projection="LL",zone=zone))
                singleCellArea<-polyArea$area
              } else {singleCellArea<-NA}
              obj@data$cellArea[x]<-singleCellArea
            }
          }
        }
        return(obj)}
         
         
