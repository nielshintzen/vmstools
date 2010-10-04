## tacsatMinimumConvexPolygon.r
## by Fabrizio Manco, 22/09/2010
## Flag vms pings inside a convex polygon regrouping a threshold percentage of total points (default 90%, for indicator DCF6)

tacsatMinimumConvexPolygon <- function (tacsat, 
                                        pctThreshold=90)
{
  require(aspace)

  vmsPoints<-cbind(tacsat$SI_LONG, tacsat$SI_LATI)

  vmsMCP<-calc_mcp(id=1, points = vmsPoints, filename="", verbose = FALSE, pct = pctThreshold)
  #plot_mcp(plotnew=TRUE, plotpoints=TRUE, titletxt="Title", xaxis= "Easting (m)", yaxis="Northing (m)")

  MCPolygon<-as.data.frame(vmsMCP$MCP.coords)
  pointInOut<-point.in.polygon(vmsPoints[,1], vmsPoints[,2], MCPolygon[,2], MCPolygon[,3])

  tacsat$INMCP<-pointInOut
  
  return(tacsat)
}