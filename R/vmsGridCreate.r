#vmsGridCreate.r
#andy south 10/2/09

#flexible function to create fishing activity grids from VMS
#change to check upload Paris june 10

vmsGridCreate <- function( dF
                         , nameLon = "Longitude"
                         , nameLat = "Latitude"
                         , nameVarToSum = ""
                         , cellsizeX = 0.5
                         , cellsizeY = 0.5
                         , we=""
                         , ea=""
                         , so=""
                         , no=""
                         , gridValName="fishing"
                         , plotMap = TRUE
                         , plotTitle = ""
                         , numCats = 5
                         , paletteCats = "heat.colors"
                         , addLegend = TRUE
                         , legendx='bottomleft'
                         , legendncol = 1
                         , legendtitle = "fishing activity"
                         , plotPoints = TRUE
                         , legPoints = FALSE
                         , colPoints = 1
                         , colLand = 'sienna'
                         , addICESgrid = FALSE
                         , addScale = TRUE
                         , outGridFile = ""  #name for output gridAscii
                         , outPlot = ""  #name for output png
                         , ... )
{

require(sp)
require(maptools)

lstargs <- list(...)

#only create grids when num records >0 (otherwise generates error)
if ( nrow(dF) > 0 )
   {
   
    #if bounds are not specified then set them from the data
    #rounds bounds to nearest whole cell unit
    if ( we == "" ) {we = min( dF[[nameLon]], na.rm=TRUE ); we = we - we%%cellsizeX}
    if ( ea == "" ) {ea = max( dF[[nameLon]], na.rm=TRUE ); ea = ea - ea%%cellsizeX + cellsizeX}
    if ( so == "" ) {so = min( dF[[nameLat]], na.rm=TRUE ); so = so - so%%cellsizeY}
    if ( no == "" ) {no = max( dF[[nameLat]], na.rm=TRUE ); no = no - no%%cellsizeY + cellsizeY}
    
    #if ( ea == "" ) ea = max( dF[[nameLon]], na.rm=T )
    #if ( so == "" ) so = min( dF[[nameLat]], na.rm=T )
    #if ( no == "" ) no = max( dF[[nameLat]], na.rm=T )

    #this copes with negative we or so values
    numXcells <- ceiling((ea-we)/cellsizeX)
    numYcells <- ceiling((no-so)/cellsizeY)
    #this copes with negative ea or no values
    numXcells <- abs(numXcells)
    numYcells <- abs(numYcells)

    #setting grid topology using package 'sp'
    #gridTopology <- GridTopology(c(we,so), c(cellsize,cellsize), c(numXcells,numYcells))
    #this sets grid at lower left corner rather than centre
    gridTopology <- GridTopology(c(we+(cellsizeX/2),so+(cellsizeY/2)), c(cellsizeX,cellsizeY), c(numXcells,numYcells))
    spatialGrid <- SpatialGrid(grid=gridTopology)
    gridded(spatialGrid) = TRUE

    #put points into a 'SpatialPointsDataFrame' 'sp' object
    coords <- cbind(x=dF[[nameLon]],y=dF[[nameLat]])
    sPDF <- SpatialPointsDataFrame(coords,data=dF)

    #overlay to find which grid cell that each VMS point is in
    #can take a long time, but does work eventually
    gridCellIndexPerVMSpoint <- over( as(sPDF,"SpatialPoints"),spatialGrid )
    sPDF$gridCellIndex <- gridCellIndexPerVMSpoint

    #if there's a column of time intervals then sum them
    if (nameVarToSum != "")
       {
        #sum timeInterval for each gridCellIndex
        perCell <- tapply(sPDF[[nameVarToSum]], gridCellIndexPerVMSpoint, sum)
        #print('1')
       } 
       else
       {
        #if no time interval, just count pings
        #probably a better way than doing this!
        sPDF$ones <- 1 #this just creates a vector of all 1s
        perCell <- tapply(sPDF$ones, gridCellIndexPerVMSpoint, sum)
        #print('2')
       }
   
    #then need to get those aggregated data values back onto the original grid
    #get it to a spatialGridDataFrame to be able to do that

    #create blank dataframe for data based on grid dimensions I know : numXcells * numYcells
    dFdataEmpty <- data.frame( seqID=(1:(numXcells*numYcells)))
    spatialGridDataFrame <- SpatialGridDataFrame( grid=gridTopology, data=dFdataEmpty )

    spatialGridDataFrame[[gridValName]] <- NA #setting up blank column
    #assigns summed values per cell back to the grid
    #!bit tricky
    spatialGridDataFrame[[gridValName]][as.integer(names(perCell))] <- c(perCell)

    if (plotMap)
       {
       
        mapGrid(spatialGridDataFrame, sPDF
               ,we=we, ea=ea, so=so, no=no
               ,gridValName=gridValName, plotTitle = plotTitle
               ,numCats = numCats, paletteCats =paletteCats, addLegend = addLegend
               ,legendx=legendx, legendncol = legendncol
               ,legendtitle = legendtitle, plotPoints = plotPoints, legPoints = legPoints, colPoints=colPoints 
               , colLand=colLand, addICESgrid = addICESgrid, addScale = addScale
               ,outGridFile = outGridFile, outPlot = outPlot, breaks0=lstargs$breaks0 )
       } #end of plotMap

    #to output spatialGridDataFrame to a gridascii file
    if (outGridFile != "")
       {
        writeAsciiGrid(spatialGridDataFrame, outGridFile, na.value = -99.99, attr=gridValName, dec='.')
       }

    #option to save the plot
    if (outPlot != "")
       {
        savePlot(outPlot,type='png')
       }
       
    #returning invisibly
    invisible(spatialGridDataFrame)   
       
   } #end of if (nrow(dF)>0)
            
} #end of vmsGridCreate



