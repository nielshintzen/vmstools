\name{indicators}
\alias{indicators}
\title{
Calculate the DCF indicators
}
\description{
This function estimates the DCF indicators 5, 6 and 7 from the tacsat dataset.
}
\usage{
indicators(indicatorNum=5,
           tacsat,
           minThreshold=10,
           pctThreshold=90,
           ltGear="",
           inShapeArea="", 
           cellresX=0.05, 
           cellresY=0.05, 
           calcAreaMethod="Trapezoid", 
           plotMapTF=FALSE,
           exportGridName="",
           exportTableName="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indicatorNum}{The indicator's number (5,6 or 7)}
  \item{tacsat}{The vms dataframe with tacsat format}
  \item{minThreshold}{The threshold value to consider a cell being "fished" in time if the ping time interval has been calculated and named SI_INTV, otherwise in number of pings}
  \item{pctThreshold}{The threshold value representing the percentage of points to include in the Minimal Convex Polygon for Indicator 6}
  \item{ltGear}{The list of gear codes to consider, if a gear code is present in the tacsat dataframe and named LE_GEAR}
  \item{inShapeArea}{The input shapefile to consider for DCF Indicator 7; path and namefile without the extension (.shp)}
  \item{cellresX}{The cell size along axis X}
  \item{cellresY}{The cell size along axis Y}
  \item{calcAreaMethod}{The method used to calculate the cell area, can be "Trapezoid" (quick and less acurate) or "UTM" (slow and accurate)}
  \item{plotMapTF}{Plot the maps}
  \item{exportGridName}{If mentionned, each grid will be exported as a ASCII grid}
  \item{exportTableName}{Name of the csv file containing the results}
}
\details{
SUMMARY

The EU Data Collection Framework (DCF) standardizes three indicators to analyse the fishing activity. They are summarised as follow:

DCF Indicator 5: Distribution of fishing activities. The spatial extent of fishing activity based on the total area of grids within which VMS records were obtained, each month;
DCF Indicator 6: Aggregation of fishing activities. The extent to which fishing activity is aggregated based on the total area of grids within which 90 percent of VMS records were obtained, each month.
DCF Indicator 7: Areas not impacted by mobile bottom gears. The area of seabed that has not been impacted by mobile bottom fishing gears in the last year. Could be reported annually and would state the total proportion of the area by depth strata in each marine region.

METHODS

These indicators aggregate the tacsat point data into a gridded data frame using the functions mapGrid.r and vmsGridCreate.r and therefore the resolution of the grid (cell size) must be defined.

DCF 5 calculates the total area of a grid of cells with fishing activity which is above a minimum threshold of number of pings or number of fishing hours (if the tacsat data contains a field with time interval between two points called SI_INTV, then the threshold will be a minimal number of hours, otherwise it will be a minimal number of points). The area of each cell is calculated with the function surface.r either via a fast and rough method using a trapezoid approximation (option "Trapezoid"), either via a more precise but slow method using a Universal Transverse Mercator projection (option "UTM"). The first method is fine for big grids of small cell sizes, the second one is better for large cell sizes. This total fishing area is processed by month.

DCF 6 also calculates the total area of a grid with fishing activity but keeps only the 90 percent of the points by discarding the outer 10\% points (or any other specified percentage). It uses the function tacsatMCP.r adapted from the aspace library. This function draws a minimum convex polygon around the central points to keep. Then these points are gridded and the total area of the cells is calculated with the surface.r function with the same optional methods as DCF 5. This total fishing area is processed by month.

DCF 7 calculates the total area of a specified polygon not impacted by mobile bottom gear. It therefore needs that the tacsat data has been merged with the logbooks in order to have a gear code (or others) for each vms point. The indicator needs a list of gear code to include as mobile bottom gears (if empty, all the points will be included). The specified area to be processed is a polygon shapefile. This polygon (or group of polygons) is then clipped with the Europe polygon to be sure that the indicator won't include land in its area calculation. If no shapefile is defined, the area of the bounding box containing all the vms points will be considered. The result is the area of the polygon less the area of the grid where fishing activity occurs. The vms pings are gridded with an optional threshold in either minimal of fishing hours or minimal number of points (see DCF 5). The area of each grid cell is calculated with the surface.r function (see DCF 5 or DCF 6).
}
\value{
For DCF Indicator 5: a list of monthly areas is returned, saved (if exportTableName is populated) and monthly grids are exported (if exportGridName is populated)
For DCF Indicator 6: a list of monthly areas is returned, saved (if exportTableName is populated) and monthly grids are exported (if exportGridName is populated)
For DCF Indicator 7: a annual value is returned and a grid is exported (if exportGridName is populated) 
}
\references{EU lot 2 project}
\author{Fabrizio Manco}
\seealso{\code{\link{mapGrids}}
         \code{\link{vmsGridCreate}}
         \code{\link{tacsatMCP}}
         \code{\link{surface}}
}
\examples{
# load the library
library(vmstools)
# Load the tacsat data
data(tacsat)

# process the tacsat data:
# mandatory if you want the gridding based on a time threshold
#   (minimal number of hours)
#pointInHarbour.r
#filterTacsat.r
#intervalTacsat.r


# load the eflalo
data(eflalo)
# merge eflalo and tacsat # mandatory for DCF Indicator 7 to
#   consider only a list of gear codes

# DCF Indicator 5
indicators(indicatorNum=5,
           tacsat,
           minThreshold=0,
           cellresX=0.05,
           cellresY=0.05,
           calcAreaMethod="Trapezoid",
           plotMapTF=TRUE,
           exportTableName="",
           exportGridName="")

# DCF Indicator 6
indicators(indicatorNum=6,
           tacsat,
           pctThreshold=90,
           cellresX=0.05,
           cellresY=0.05,
           calcAreaMethod="Trapezoid",
           plotMapTF=TRUE,
           exportTableName="",
           exportGridName="")

# DCF Indicator 7
\dontrun{
indicators(indicatorNum=7,
           tacsat, 
           ltGear=c("TBB","OTB","PTB","DRB","DRH"),
           inShapeArea="Shapefile",
           cellresX=0.05,
           cellresY=0.05, 
           calcAreaMethod="Trapezoid",
           minThreshold=0,
           plotMapTF=TRUE,
           exportGridName="")
}
}