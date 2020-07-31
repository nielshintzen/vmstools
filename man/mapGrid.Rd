\name{mapGrid}
\alias{mapGrid}
\title{function to map grids}
\description{
Accepts an input of a \code{SpatialGridDataFrame}
Plots a map of the grid and optionally outputs to a gridAscii file and/or an image. 
}
\usage{
mapGrid(sGDF
             , sPDF
             , we = ""
             , ea = ""
             , so = ""
             , no = ""
             , gridValName="fishing"
             , plotTitle = ""
             , numCats = 5
             , paletteCats = "heat.colors"
             , addLegend = TRUE
             , legendx = "bottomleft"
             , legendncol = 1
             , legendtitle = "fishing activity"
             , plotPoints = FALSE
             , colPoints = 1
             , legPoints = FALSE
             , colland = "sienna"
             , addICESgrid = FALSE
             , addScale = TRUE
             , outGridFile = ""  
             , outPlot = ""
             , ...)
}

\arguments{
  \item{sGDF}{a \code{SpatialGridDataFrame}}
  \item{sPDF}{an optional \code{SpatialPointsDataFrame} plotted if \code{plotPoints=TRUE}}
  \item{we}{western bounds of the area to plot, if not specified taken from the \code{sGDF}}
  \item{ea}{eastern bounds of the area to plot, if not specified taken from the \code{sGDF}}
  \item{so}{southern bounds of the area to plot, if not specified taken from the \code{sGDF}}
  \item{no}{northern bounds of the area to plot, if not specified taken from the \code{sGDF}}
  \item{gridValName}{the name of the attribute column to plot from the \code{SpatialGridDataFrame}}
  \item{plotTitle}{optional title to add to the plot}
  \item{numCats}{how many categories to classify grid values into for map plot (uses\code{pretty()}) classification)}
  \item{paletteCats}{color pallete to use}
  \item{addLegend}{whether to add a legend to the plot}
  \item{legendx}{position of legend should be one of 'bottomright', 'bottom', 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right', 'center'}
  \item{legendncol}{number of columns in the legend}
  \item{legendtitle}{legend title}
  \item{plotPoints}{whether to add the original points to the plot}
  \item{colPoints}{color of points to plot}
  \item{legPoints}{Logical. Points in legend}
  \item{colland}{color of land}
  \item{addICESgrid}{Logical. Adding ICES grid on top}
  \item{addScale}{Logical. Adding axes}
  \item{outGridFile}{optional name for a gridAscii file to be created from the grid}
  \item{outPlot}{optional name for a png file to be created from the plot}
  \item{\dots}{NOT used yet}
}
\references{EU Lot 2 project}
\author{Andy South}
\seealso{\code{vmsGridCreate()}}
\examples{
#mapGrid(dF, nameLon = "POS_LONGITUDE", nameLat = "POS_LATITUDE",
#        cellsizeX = 0.5, cellsizeY = 0.5,legendx='bottomright',
#        plotPoints=TRUE )
}

