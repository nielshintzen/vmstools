#vmsGridCreate.r
#andy south 10/2/09

#flexible function to create fishing activity grids from VMS
#change to check upload Paris june 10



#' function to create grids from point data by counting points in cells or
#' summing an attribute
#' 
#' Accepts an input of points data in a data frame with named columns for X and
#' Y. Creates a grid of defined cell size in X and Y directions (cells can be
#' unequal). Either counts points in cells or summs a attribute variable if one
#' is supplied. Optionally plots a map of the grid and outputs to a gridAscii
#' file and/or an image.
#' 
#' 
#' @param dF a dataFrame containing point data
#' @param nameLon name of the column in the dataFrame containing Longitude or x
#' values
#' @param nameLat name of the column in the dataFrame containing Latitude or y
#' values
#' @param nameVarToSum optional name of the column in the dataFrame containing
#' the attribute values to sum in the grid. If set to "" points are counted
#' @param cellsizeX length X (horizontal) of desired grid cells, in same units
#' as the coordinates
#' @param cellsizeY length Y (vertical) of desired grid cells, in same units as
#' the coordinates
#' @param we western bounds of the desired grid
#' @param ea eastern bounds of the desired grid
#' @param so southern bounds of the desired grid
#' @param no northern bounds of the desired grid
#' @param gridValName the name to give to the attribute column of the returned
#' \code{SpatialGridDataFrame}, set to 'fishing' by default
#' @param plotMap whether to plot a map of the resulting grid
#' @param plotTitle optional title to add to the plot
#' @param numCats how many categories to classify grid values into for map plot
#' (uses\code{pretty()}) classification)
#' @param paletteCats color pallete to use
#' @param addLegend whether to add a legend to the plot
#' @param legendx position of legend should be one of 'bottomright', 'bottom',
#' 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right', 'center'
#' @param legendncol number of columns in the legend
#' @param legendtitle legend title
#' @param plotPoints whether to add the original points to the plot
#' @param legPoints Logical. Points in legend
#' @param colPoints color of points to plot
#' @param colland color of land
#' @param addICESgrid Logical. Adding ICES grid on top
#' @param addScale Logical. Adding axes
#' @param outGridFile optional name for a gridAscii file to be created from the
#' grid
#' @param outPlot optional name for a png file to be created from the plot
#' @param \dots NOT used yet
#' @return a \code{SpatialGridDataFrame} object of the grid defined in package
#' \code{sp}
#' @author Andy South
#' @seealso \code{\link{mapGrid}}
#' @references EU VMS tools project
#' @examples
#' 
#' #vmsGridCreate(dF, nameLon = "POS_LONGITUDE", nameLat = "POS_LATITUDE",
#' # cellsizeX = 0.5, cellsizeY = 0.5,legendx='bottomright',plotPoints=TRUE )
#' #get the example data
#' data(tacsat)
#' 
#' #subset the first 2000 points to avoid problems with NAs
#' dFVMS <- tacsat[1:2000,]
#' 
#' #create vms grid minimum call with defaults
#' vmsGridCreate(dFVMS,nameLat='SI_LATI',nameLon='SI_LONG')
#' 
#' #making the grid finer
#' vmsGridCreate(dFVMS,nameLat='SI_LATI',nameLon='SI_LONG',
#'               cellsizeX=0.05,cellsizeY=0.05)
#' 
#' 
#' @export vmsGridCreate
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



