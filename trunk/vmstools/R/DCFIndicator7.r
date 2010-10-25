## DCFIndicator7.r
## by Fabrizio Manco, 23/09/2010
## calculates the DCF7 indicator : area of cells not containing activity

DCFIndicator7 <- function ( tacsat,
                            mobileBottomGear="",          # a list of gear code
                            inShapeArea="",               # the name of the shapefile without the .shp extension
                            cellresX=0.05,
                            cellresY=0.05,
                            calcAreaMethod="Trapezoid",   # "Trapezoid" (fast and less accurate, good for small cellsizes) or "UTM" (accurate but slow, good for huge cellsizes)
                            minThreshold=10,              # if time interval has been calculated (and named SI_INTV), it's a minimal nb of minutes, otherwise, it's minimal number of points
                            plotMapTF = FALSE
                            )          
                                                 
{ require(shapefiles)
  require(sp)
  require(PBSmapping)
  tacsat<-tacsat[complete.cases(tacsat),]
  if (mobileBottomGear!="") {tacsat<-subset(tacsat, tacsat$LE_GEAR %in% mobileBottomGear)}

  if (inShapeArea!="")
    { # read the shapefile 

      shapeAll<-read.shapefile(inShapeArea)
    
      # clip the shape polygon with the land
      clipShapeFromLand<-clipPolygons (shapeAll, europa)
    
      vmsPingsCoord<-cbind(tacsat$SI_LONG, tacsat$SI_LATI)
      pointInOutByPoly<-rep(0,length(vmsPingsCoord[,1]))
      
      ltPoly<-unique(clipShapeFromLand$PID)

      # points in polygons
      for (x in 1:length(ltPoly)){                                   
        polyCoord<-cbind(clipShapeFromLand$X[clipShapeFromLand$PID==ltPoly[x]],clipShapeFromLand$Y[clipShapeFromLand$PID==ltPoly[x]])
        pointInOutByPoly<-pointInOutByPoly + point.in.polygon(vmsPingsCoord[,1], vmsPingsCoord[,2], polyCoord[,1], polyCoord[,2])
        }

      tacsat$pointInOut<-pointInOutByPoly
      tacsat<-subset(tacsat, pointInOut!=0)
      }
  
  # Grid the points
  if ("SI_INTV" %in% colnames(tacsat)) { nameVarToSum="SI_INTV"} else {nameVarToSum=""}
  vmsGrid<-vmsGridCreate(tacsat, nameLon = "SI_LONG", nameLat = "SI_LATI", cellsizeX=cellresX, cellsizeY=cellresY, nameVarToSum, plotMap=plotMapTF, plotPoints = FALSE)
  
  # calculate the area of each cell in square km
  vmsGrid<-surface(vmsGrid, method=calcAreaMethod, includeNA=TRUE)
                                                                     
  if (inShapeArea!="")
    {    
      # specify which grid cell is in the polygon
      gridPointsCoord<-coordinates(vmsGrid)
      gridCellInOutByPoly<-rep(0,length(gridPointsCoord[,1]))
      
      # cells in polygon
      for (x in 1:length(ltPoly)){
        polyCoord<-cbind(clipShapeFromLand$X[clipShapeFromLand$PID==ltPoly[x]],clipShapeFromLand$Y[clipShapeFromLand$PID==ltPoly[x]])
        gridCellInOutByPoly<-gridCellInOutByPoly + point.in.polygon(gridPointsCoord[,1], gridPointsCoord[,2], polyCoord[,1], polyCoord[,2])

        }
        
      vmsGrid$inPolygon<-gridCellInOutByPoly
      areaInPolygon<-sum(vmsGrid@data$cellArea[vmsGrid@data$inPolygon==1])
      
      } else {areaInPolygon<-sum(vmsGrid@data$cellArea)}
      
  # calculate the areas  
  areaFishing<-sum(vmsGrid@data$cellArea[!is.na(vmsGrid@data$fishing) & vmsGrid@data$fishing>minThreshold])
  
  # calculate the result
  resultDCF7<-areaInPolygon-areaFishing
  
  return(resultDCF7)
}  