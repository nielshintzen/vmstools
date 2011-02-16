\name{mosaicASCIIGrids}
\alias{mosaicASCIIGrids}
\title{
Merged a list of grids
}
\description{
This function merges a list of ASCII grid files with different extents together and export it as a new grid.
}
\usage{
mosaicASCIIGrids(ascgridfilelist, outputFileName)
}
\arguments{
  \item{ascgridfilelist}{a list of ASCII grids filenames}
  \item{outputFileName}{the name of the resulting merged grid}
}
\details{
This function is based on the mosaic function of the raster package. It will merge a serie of grids together, which might have been created by vmsGridCreate. 
The values of overlapping cells is summed.
}
\value{
The result is a new ASCII grid file.
}
\references{EU lot 2 project}
\author{Fabrizio Manco}
\seealso{\code{indicators()}}
\examples{
# Each country can create a ASCII grid using vmsGridCreate, for example
# Uk
UkVmsGrid<-vmsGridCreate(UkTacsat, nameLon = "SI_LONG", nameLat = "SI_LATI", cellsizeX=0.05, cellsizeY=0.05, nameVarToSum="", outGridFile="UkVmsGrid.asc")
# Dutch
NLVmsGrid<-vmsGridCreate(NLTacsat, nameLon = "SI_LONG", nameLat = "SI_LATI", cellsizeX=0.05, cellsizeY=0.05, nameVarToSum="", outGridFile="NLVmsGrid.asc")
# and so on with other countries...

# List of grids to merge
ascgridfilelist<-c("UkVmsGrid.asc","NLVmsGrid.asc")

mosaicASCIIGrids(ascgridfilelist, "Uk+NLVmsGrids.asc")
}
