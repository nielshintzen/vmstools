\name{vmsGridCreate}
\alias{vmsGridCreate}
\title{function to create grids from point data by counting points in cells or summing an attribute}
\description{
Accepts an input of points data in a data frame with named columns for X and Y.
Creates a grid of defined cell size in X and Y directions (cells can be unequal).
Either counts points in cells or summs a attribute variable if one is supplied.
Optionally plots a map of the grid and outputs to a gridAscii file and/or an image. 
}
\usage{
vmsGridCreate(dF
             , nameLon = "Longitude"
             , nameLat = "Latitude"
             , nameVarToSum = ""
             , cellsizeX = 0.5
             , cellsizeY = 0.5
             , we = ""
             , ea = ""
             , so = ""
             , no = ""
             , gridValName="fishing"
             , plotMap = TRUE
             , plotTitle = ""
             , numCats = 5
             , addLegend = TRUE
             , legendx = "bottomleft"
             , legendncol = 1
             , legendtitle = "fishing activity"
             , plotPoints = TRUE
             , outGridFile = ""  
             , outPlot = ""
             , ...)
}

\arguments{
  \item{dF}{a dataFrame containing point data}
  \item{nameLon}{name of the column in the dataFrame containing Longitude or x values}
  \item{nameLat}{name of the column in the dataFrame containing Latitude or y values}
  \item{nameVarToSum}{optional name of the column in the dataFrame containing the attribute values to sum in the grid. If set to "" points are counted }
  \item{cellsizeX}{length X (horizontal) of desired grid cells, in same units as the coordinates}
  \item{cellsizeY}{length Y (vertical) of desired grid cells, in same units as the coordinates}
  \item{we}{western bounds of the desired grid}
  \item{ea}{eastern bounds of the desired grid}
  \item{so}{southern bounds of the desired grid}
  \item{no}{northern bounds of the desired grid}
  \item{gridValName}{the name to give to the attribute column of the returned \code{SpatialGridDataFrame}, set to 'fishing' by default}
  \item{plotMap}{whether to plot a map of the resulting grid}
  \item{plotTitle}{optional title to add to the plot}
  \item{numCats}{how many categories to classify grid values into for map plot (uses\code{pretty()}) classification)}
  \item{addLegend}{whether to add a legend to the plot}
  \item{legendx}{position of legend should be one of 'bottomright', 'bottom', 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right', 'center'}
  \item{legendncol}{number of columns in the legend}
  \item{legendtitle}{legend title}
  \item{plotPoints}{whether to add the original points to the plot}
  \item{outGridFile}{optional name for a gridAscii file to be created from the grid}
  \item{outPlot}{optional name for a png file to be created from the plot}
  \item{\dots}{NOT used yet}
}
\details{}
\value{ a \code{SpatialGridDataFrame} object of the grid defined in package \code{sp} }
\references{EU VMS tools project}
\author{Andy South}
\note{}
\seealso{\code{mapGrid()}}
\examples{
#vmsGridCreate(dF, nameLon = "POS_LONGITUDE", nameLat = "POS_LATITUDE",
# cellsizeX = 0.5, cellsizeY = 0.5,legendx='bottomright',plotPoints=T )
#get the example data
data(tacsat)
#subset the first 2000 points to avoid problems with NAs
dFVMS <- tacsat[1:2000,]
#create vms grid minimum call with defaults
vmsGridCreate(dF2,nameLat='declat',nameLon='declon')
#making the grid finer
vmsGridCreate(dF2,nameLat='declat',nameLon='declon',
              cellsizeX=0.05,cellsizeY=0.05)
#summing a data column and adding a title
vmsGridCreate(dF2,nameLat='declat',nameLon='declon',
              cellsizeX=0.05,cellsizeY=0.05,nameVarToSum='heading',
              plotTitle='lot2 paris')
#changing the number of categories in the plot
vmsGridCreate(dF2,nameLat='declat',nameLon='declon',
              cellsizeX=0.05,cellsizeY=0.05,numCats=10)
#outputting to a gridAscii file that can be loaded into ArcGIS
vmsGridCreate(dF2,nameLat='declat',nameLon='declon',
              cellsizeX=0.05,cellsizeY=0.05,outGridFile='test.asc')
}
