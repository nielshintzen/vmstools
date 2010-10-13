## calcAreaOfCells.r
## by Fabrizio Manco, 15/09/2010
## calculates the area in square kilometre of each cell of a grid
## based on PBSmapping library, each cell is projected with a UTM projection
## very accurate but very slow method

calcAreaOfCells <- function ( vmsGrid,
                              includeNA=FALSE )
{
  require(PBSmapping)

  if (class(vmsGrid)=='SpatialGridDataFrame') # not empty...
  {
    griddims <- summary(vmsGrid)$grid
    sizelon  <- griddims[1,2]
    sizelat  <- griddims[2,2]

    ltCentreCell<-coordinates(vmsGrid)

    for (x in 1:(length(ltCentreCell)/2))
      { 
      if (!is.na(vmsGrid@data[x,2]) | includeNA) {        # speed up the calculation by dropping cells with fishing=NA  /!| only work for DCF5 and DCF6!
        minX<-ltCentreCell[x,1]-sizelon/2
        maxX<-ltCentreCell[x,1]+sizelon/2
        minY<-ltCentreCell[x,2]-sizelat/2
        maxY<-ltCentreCell[x,2]+sizelat/2
        ltX<-c(minX,minX,maxX,maxX)
        ltY<-c(minY,maxY,maxY,minY)
        sourcePoly<-cbind(rep(1,4),seq(1,4),ltX,ltY)
        rownames(sourcePoly)<-seq(1,4)
        colnames(sourcePoly)<-c("PID","POS","X","Y")

        polyArea<-calcArea(as.PolySet(sourcePoly, projection="LL"))
        singleCellArea<-polyArea$area
        } else {singleCellArea<-NA}
        vmsGrid@data$cellArea[x]<-singleCellArea
      }
  }
  return(vmsGrid)
}                                                                 