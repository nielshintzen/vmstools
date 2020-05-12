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
\seealso{\code{\link{indicators}}}
\examples{
data(tacsat)
tacsat <- sortTacsat(tacsat)
# Each country can create a ASCII grid using vmsGridCreate, for example
# year 1800 from tacsat file
VmsGrid00<-vmsGridCreate(subset(tacsat,year(SI_DATIM)==1800), nameLon = "SI_LONG",
            nameLat = "SI_LATI",cellsizeX=0.05, cellsizeY=0.05,
            nameVarToSum="",outGridFile="VmsGrid00.asc")
# year 1801 from tacsat file
VmsGrid01<-vmsGridCreate(subset(tacsat,year(SI_DATIM)==1800), nameLon = "SI_LONG",
             nameLat = "SI_LATI",cellsizeX=0.05, cellsizeY=0.05,
             nameVarToSum="",outGridFile="VmsGrid01.asc")

# and so on with other countries or years...

# List of grids to merge
ascgridfilelist<-c("VmsGrid00.asc","VmsGrid01.asc")

mosaicASCIIGrids(ascgridfilelist, "VmsGrids.asc")
}
