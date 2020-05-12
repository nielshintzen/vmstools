## clipPolygons.r
## by Fabrizio Manco, 13/10/2010
## clip two polygons, used to remove lands some study polygon area for DCF7

clipPolygons <- function ( shapeAll,
                           europa)
  {
  require(adehabitat)
  require(PBSmapping)
  
  if (!exists("europa")) {load("europa.rda")}

  allSourcePoly<-numeric()

  # transform the shape data to a format compatible with PolySets
  for (x in 1:length(shapeAll$shp$shp))
    {
    eachPolyFromShape<-shapeAll$shp$shp[[x]]
    sourcePoly<-cbind(eachPolyFromShape$record,seq(1:length(eachPolyFromShape$points$X)),eachPolyFromShape$points)
    allSourcePoly<-rbind(allSourcePoly,sourcePoly)
    }

  colnames(allSourcePoly)<-c("PID","POS","X","Y")
  allPolyFromShape<-as.PolySet(allSourcePoly)
  
  # clip the polygons
  clipShapeFromLand<-joinPolys(allPolyFromShape,europa, operation="DIFF",maxVert=1e+06)
  
  # rearrange the polygons to remove the sub-polygons
  clipShapeFromLand<-combinePolys(clipShapeFromLand)
  clipShapeFromLand<-cbind(clipShapeFromLand$SID, clipShapeFromLand$POS, clipShapeFromLand$X, clipShapeFromLand$Y)
  colnames(clipShapeFromLand)<-c("PID","POS","X","Y")
  clipShapeFromLand<-as.PolySet(clipShapeFromLand)
  
  return(clipShapeFromLand)
  }



