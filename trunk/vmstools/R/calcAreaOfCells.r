## calcAreaOfCells.r
## by Fabrizio Manco, 15/09/2010
## calculates the area in square kilometre of each cell of a grid

calcAreaOfCells <- function (vmsGrid)
{
  #require(degree2Km.r)

  if (class(vmsGrid)=='SpatialGridDataFrame') # not empty...
  {
    gridCoord<-coordinates(vmsGrid)
    gridpar<-gridparameters(vmsGrid)
    cellSizeX<-gridpar$cellsize[1]
    cellSizeY<-gridpar$cellsize[2]

    gridCoordMinX<-gridCoord[,1]-(cellSizeX/2)
    gridCoordMaxX<-gridCoord[,1]+(cellSizeX/2)
    gridCoordMinY<-gridCoord[,2]-(cellSizeY/2)
    gridCoordMaxY<-gridCoord[,2]+(cellSizeY/2)

    # method 1 using the degree2Km script
    gridSizeKmY<-(abs(gridCoordMaxY-gridCoordMinY))*111.12
    gridSizeKmX<-degree2Km(gridCoordMinX,gridCoordMinY,cellSizeX)
    
    vmsGrid@data$cellArea<-gridSizeKmX*gridSizeKmY
    
    return(vmsGrid)
    }
}