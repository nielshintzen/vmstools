## tacsatMinimumConvexPolygon.r
## by Fabrizio Manco, 22/09/2010
## Flag vms pings inside a convex polygon regrouping a threshold percentage of total points (default 90%, for indicator DCF6)



#' Flag tacsat records that are within the convex polygon
#' 
#' Flag tacsat records that are within the convex polygon with a predefined
#' threshold
#' 
#' See point.in.polygon function for returned value details
#' 
#' @param tacsat Tacsat dataframe
#' @param pctThreshold Threshold of points to consider. Between 0 and 100.
#' @author Fabrizio Manco
#' @references EU Lot 2 project
#' @examples
#' 
#' require(adehabitat)
#' require(ade4)
#' data(tacsat)
#' tacsat <- tacsat[1:100,]
#' tacsatMinimumConvexPolygon(tacsat,95)
#' 
#' @export tacsatMinimumConvexPolygon
tacsatMinimumConvexPolygon <- function (tacsat, 
                                        pctThreshold=90)
{
  require(aspace)
  require(ade4)

  vmsPoints<-cbind(tacsat$SI_LONG, tacsat$SI_LATI)

  vmsMCP<-calc_mcp(id=1, points = vmsPoints, filename="", verbose = FALSE, pct = pctThreshold)

  MCPolygon<-as.data.frame(vmsMCP$MCP.coords)
  pointInOut<-point.in.polygon(vmsPoints[,1], vmsPoints[,2], MCPolygon[,2], MCPolygon[,3])

  tacsat$INMCP<-pointInOut
  
  return(tacsat)
}
