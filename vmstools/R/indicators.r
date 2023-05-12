## indicators.r
## by Fabrizio Manco, 14/02/2011
## calculates the DCF indicators 5,6 or 7



#' Calculate the DCF indicators
#' 
#' This function estimates the DCF indicators 5, 6 and 7 from the tacsat
#' dataset.
#' 
#' SUMMARY
#' 
#' The EU Data Collection Framework (DCF) standardizes three indicators to
#' analyse the fishing activity. They are summarised as follow:
#' 
#' DCF Indicator 5: Distribution of fishing activities. The spatial extent of
#' fishing activity based on the total area of grids within which VMS records
#' were obtained, each month; DCF Indicator 6: Aggregation of fishing
#' activities. The extent to which fishing activity is aggregated based on the
#' total area of grids within which 90 percent of VMS records were obtained,
#' each month. DCF Indicator 7: Areas not impacted by mobile bottom gears. The
#' area of seabed that has not been impacted by mobile bottom fishing gears in
#' the last year. Could be reported annually and would state the total
#' proportion of the area by depth strata in each marine region.
#' 
#' METHODS
#' 
#' These indicators aggregate the tacsat point data into a gridded data frame
#' using the functions mapGrid.r and vmsGridCreate.r and therefore the
#' resolution of the grid (cell size) must be defined.
#' 
#' DCF 5 calculates the total area of a grid of cells with fishing activity
#' which is above a minimum threshold of number of pings or number of fishing
#' hours (if the tacsat data contains a field with time interval between two
#' points called SI_INTV, then the threshold will be a minimal number of hours,
#' otherwise it will be a minimal number of points). The area of each cell is
#' calculated with the function surface.r either via a fast and rough method
#' using a trapezoid approximation (option "Trapezoid"), either via a more
#' precise but slow method using a Universal Transverse Mercator projection
#' (option "UTM"). The first method is fine for big grids of small cell sizes,
#' the second one is better for large cell sizes. This total fishing area is
#' processed by month.
#' 
#' DCF 6 also calculates the total area of a grid with fishing activity but
#' keeps only the 90 percent of the points by discarding the outer 10\% points
#' (or any other specified percentage). It uses the function tacsatMCP.r
#' adapted from the aspace library. This function draws a minimum convex
#' polygon around the central points to keep. Then these points are gridded and
#' the total area of the cells is calculated with the surface.r function with
#' the same optional methods as DCF 5. This total fishing area is processed by
#' month.
#' 
#' DCF 7 calculates the total area of a specified polygon not impacted by
#' mobile bottom gear. It therefore needs that the tacsat data has been merged
#' with the logbooks in order to have a gear code (or others) for each vms
#' point. The indicator needs a list of gear code to include as mobile bottom
#' gears (if empty, all the points will be included). The specified area to be
#' processed is a polygon shapefile. This polygon (or group of polygons) is
#' then clipped with the Europe polygon to be sure that the indicator won't
#' include land in its area calculation. If no shapefile is defined, the area
#' of the bounding box containing all the vms points will be considered. The
#' result is the area of the polygon less the area of the grid where fishing
#' activity occurs. The vms pings are gridded with an optional threshold in
#' either minimal of fishing hours or minimal number of points (see DCF 5). The
#' area of each grid cell is calculated with the surface.r function (see DCF 5
#' or DCF 6).
#' 
#' @param indicatorNum The indicator's number (5,6 or 7)
#' @param tacsat The vms dataframe with tacsat format
#' @param minThreshold The threshold value to consider a cell being "fished" in
#' time if the ping time interval has been calculated and named SI_INTV,
#' otherwise in number of pings
#' @param pctThreshold The threshold value representing the percentage of
#' points to include in the Minimal Convex Polygon for Indicator 6
#' @param ltGear The list of gear codes to consider, if a gear code is present
#' in the tacsat dataframe and named LE_GEAR
#' @param inShapeArea The input shapefile to consider for DCF Indicator 7; path
#' and namefile without the extension (.shp)
#' @param cellresX The cell size along axis X
#' @param cellresY The cell size along axis Y
#' @param calcAreaMethod The method used to calculate the cell area, can be
#' "Trapezoid" (quick and less acurate) or "UTM" (slow and accurate)
#' @param plotMapTF Plot the maps
#' @param exportGridName If mentionned, each grid will be exported as a ASCII
#' grid
#' @param exportTableName Name of the csv file containing the results
#' @return For DCF Indicator 5: a list of monthly areas is returned, saved (if
#' exportTableName is populated) and monthly grids are exported (if
#' exportGridName is populated) For DCF Indicator 6: a list of monthly areas is
#' returned, saved (if exportTableName is populated) and monthly grids are
#' exported (if exportGridName is populated) For DCF Indicator 7: a annual
#' value is returned and a grid is exported (if exportGridName is populated)
#' @author Fabrizio Manco
#' @seealso \code{\link{mapGrids}} \code{\link{vmsGridCreate}}
#' \code{\link{tacsatMCP}} \code{\link{surface}}
#' @references EU lot 2 project
#' @examples
#' 
#' # load the library
#' library(vmstools)
#' # Load the tacsat data
#' data(tacsat)
#' 
#' # process the tacsat data:
#' # mandatory if you want the gridding based on a time threshold
#' #   (minimal number of hours)
#' #pointInHarbour.r
#' #filterTacsat.r
#' #intervalTacsat.r
#' 
#' 
#' # load the eflalo
#' data(eflalo)
#' # merge eflalo and tacsat # mandatory for DCF Indicator 7 to
#' #   consider only a list of gear codes
#' 
#' # DCF Indicator 5
#' indicators(indicatorNum=5,
#'            tacsat,
#'            minThreshold=0,
#'            cellresX=0.05,
#'            cellresY=0.05,
#'            calcAreaMethod="Trapezoid",
#'            plotMapTF=TRUE,
#'            exportTableName="",
#'            exportGridName="")
#' 
#' # DCF Indicator 6
#' indicators(indicatorNum=6,
#'            tacsat,
#'            pctThreshold=90,
#'            cellresX=0.05,
#'            cellresY=0.05,
#'            calcAreaMethod="Trapezoid",
#'            plotMapTF=TRUE,
#'            exportTableName="",
#'            exportGridName="")
#' 
#' # DCF Indicator 7
#' \dontrun{
#' indicators(indicatorNum=7,
#'            tacsat, 
#'            ltGear=c("TBB","OTB","PTB","DRB","DRH"),
#'            inShapeArea="Shapefile",
#'            cellresX=0.05,
#'            cellresY=0.05, 
#'            calcAreaMethod="Trapezoid",
#'            minThreshold=0,
#'            plotMapTF=TRUE,
#'            exportGridName="")
#' }
#' 
#' @export indicators
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
