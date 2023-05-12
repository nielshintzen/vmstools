## tacsatMinimumConvexPolygon.r
## by Fabrizio Manco, 22/09/2010
## Flag vms pings inside a convex polygon regrouping a threshold percentage of total points (default 90%, for indicator DCF6)



#' Flag the tacsat ping inside a Minimum Convex Polygon
#' 
#' This function will flag the tacsat pings that are located inside a minimum
#' convex polygon which contains a percentage of the whole number of pings
#' 
#' This function is based on the calc_mcp function of the aspace package. The
#' plot_mcp function allows to plot the polygon. It is used in the DCF
#' Indicator 6 as a filtering of aggregated fisheries.
#' 
#' @param tacsat tacsat data frame
#' @param pctThreshold Percentage of pings to be included into the MCP
#' @return The tacsat is returned with attached a column called INMCP which
#' flag the pings inside the MCP. The nodes of the MCP are also returned.
#' @author Fabrizio Manco
#' @seealso \code{\link{indicators}}
#' @references EU lot 2 project
#' @examples
#' 
#' require(PBSmapping)
#' data(tacsat)
#' 
#' # Flag the pings inside a polygon gathering 90% of the pings
#' tacsat<-tacsatMCP(tacsat, pctThreshold=90)
#' # Filter the tacsat data to remove the points outside the MCP
#' tacsat<-subset(tacsat, tacsat$INMCP!=0)
#' # Grid the filtered vms points
#' vmsgrid<-vmsGridCreate(tacsat, nameLon = "SI_LONG", nameLat = "SI_LATI",
#'                        cellsizeX=0.05, cellsizeY=0.05, plotMap=TRUE)
#' # Add the MCP to the map
#' plot_mcp(plotnew=FALSE, plotpoints=FALSE, titletxt="")
#' 
#' @export tacsatMCP
tacsatMCP <- function (tacsat, pctThreshold=90)
{
  require(aspace)

  vmsPoints<-cbind(tacsat$SI_LONG, tacsat$SI_LATI)

  vmsMCP<-calc_mcp(id=1, points = vmsPoints, filename="", verbose = FALSE, pct = pctThreshold)

  MCPolygon<-as.data.frame(vmsMCP$MCP.coords)
  pointInOut<-point.in.polygon(vmsPoints[,1], vmsPoints[,2], MCPolygon[,2], MCPolygon[,3])

  tacsat$INMCP<-pointInOut
  
  return(tacsat)
}
