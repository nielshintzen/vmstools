## indicators.r
## by Fabrizio Manco, 14/02/2011
## calculates the DCF indicators 5,6 or 7

indicators <- function    ( indicatorNum=5,             # indicator 5, 6 or 7
                            tacsat,                     # tacsat-like input data
                            minThreshold=10,            # if time interval has been calculated (and named SI_INTV), it's a minimal nb of minutes, otherwise, it's minimal number of points
                            pctThreshold=90,            # specific to indicator 6, percentage of points to include
                            ltGear="",                  # a list of gear code to keep for the analysis /!\ gear code field must be called LE_GEAR
                            inShapeArea="",             # specific to indicator 7, the name of the shapefile without the .shp extension
                            cellresX=0.05,              # grid cell resolution, x axis
                            cellresY=0.05,              # grid cell resolution, y axis
                            calcAreaMethod="Trapezoid", # "Trapezoid" (fast and less accurate, good for small cellsizes) or "UTM" (accurate but slow, good for huge cellsizes)
                            plotMapTF=FALSE,
                            exportGridName="",          # if not empty, the monthly (DCF5 and 6) grids will be exported as an ASCII grid (.asc) named with this string, DCF number and month number
                            exportTableName=""
                            )
  {
  # for all indicators
  #remove rows with NA
  tacsat<-tacsat[complete.cases(tacsat),]
  # gear filtering: intially only necessary for DCF 7 ("area impacted by mobile bottom gears"), but might be useful for other indicators
  # keep only the gear codes listed in mobileBottomGear, and tacsat must have a LE_GEAR column
  if (length(ltGear)>1 & !is.null(tacsat$LE_GEAR)) {tacsat<-subset(tacsat, tacsat$LE_GEAR %in% ltGear)}
  
  if (indicatorNum==7)
    {
    #### DCF INDICATOR 7 -  Areas not impacted by mobile bottom gears ####
    
    require(shapefiles)
    require(sp)
    require(PBSmapping)

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
    if (exportGridName!="") {outGridFileName<-paste(exportGridName,"_DCF",indicatorNum,".asc",sep="")} else {outGridFileName<-""}
    vmsGrid<-vmsGridCreate(tacsat, nameLon = "SI_LONG", nameLat = "SI_LATI", cellsizeX=cellresX, cellsizeY=cellresY, nameVarToSum, plotMap=plotMapTF, plotPoints = FALSE, outGridFile=outGridFileName)

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

      } else {areaInPolygon<-sum(vmsGrid@data$cellArea)}  # if no shapefile area is provided the entire rectangle containing the pings is processed

    # calculate the areas
    areaFishing<-sum(vmsGrid@data$cellArea[!is.na(vmsGrid@data$fishing) & vmsGrid@data$fishing>minThreshold])

    # calculate the result
    tableResultDCF<-areaInPolygon-areaFishing
    }

  #### DCF INDICATORS 5 AND 6 ####
  if (indicatorNum==5 | indicatorNum==6)
    {
    # commons for DCF 5 and 6
    if (!"SI_DATIM" %in% colnames(tacsat)) {tacsat$SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M")}
    ltMonth<-unique(as.numeric(format(tacsat$SI_DATIM, format="%m")))
    ltMonth<-sort(ltMonth)
    tableResultDCF=data.frame(month=ltMonth, DCF=rep(0, length(ltMonth)))

    # monthly process
    for (x in 1:length(ltMonth)){
      currMonth<-ltMonth[x]
      monthlyTacsat<-subset(tacsat, as.numeric(format(tacsat$SI_DATIM, format="%m"))==currMonth)
      
      if (indicatorNum==5)
      {
        # specific to DCF 5
        if ("SI_INTV" %in% colnames(tacsat)) {nameVarToSum="SI_INTV"} else {nameVarToSum=""}
      }
      
      if (indicatorNum==6)
      {
        # specific to DCF 6
        # flag the pings inside the MCP according to the pctThreshold
        monthlyTacsat<-tacsatMCP(monthlyTacsat, pctThreshold)
        # grid the vms pings inside the MCP
        monthlyTacsat<-subset(monthlyTacsat, monthlyTacsat$INMCP!=0)
        nameVarToSum=""
      }
      
      # Create the grid
      if (plotMapTF) {x11()}
      if (exportGridName!="") {outGridFileName<-paste(exportGridName,"_DCF",indicatorNum,"_Month", currMonth,".asc",sep="")} else {outGridFileName<-""}
      monthlyVmsGrid<-vmsGridCreate(monthlyTacsat, nameLon = "SI_LONG", nameLat = "SI_LATI", cellsizeX=cellresX, cellsizeY=cellresY, nameVarToSum, plotMap=plotMapTF, plotTitle=paste("Month ", currMonth), plotPoints = FALSE, outGridFile=outGridFileName)
      if (plotMapTF==TRUE & indicatorNum==6) {plot_mcp(plotnew=FALSE, plotpoints=FALSE, titletxt="")}   # plot the specific DCF 6 MCP

      # calculate the area of each cell in square km and store it into a table
      monthlyVmsGrid<-surface(monthlyVmsGrid, method=calcAreaMethod, includeNA=FALSE)
      tableResultDCF[x,2]<-sum(monthlyVmsGrid@data$cellArea[!is.na(monthlyVmsGrid@data$fishing) & monthlyVmsGrid@data$fishing>minThreshold])
      }
    if (exportTableName!="") 
      {
      # export the table as a csv in the wd
      write.csv(tableResultDCF, paste(exportTableName, indicatorNum, ".csv", sep=""))
      }
    }
  return(tableResultDCF)
  }