\name{tacsatMCP}
\alias{tacsatMCP}
\title{
Flag the tacsat ping inside a Minimum Convex Polygon
}
\description{
This function will flag the tacsat pings that are located inside a minimum convex polygon which contains a percentage of the whole number of pings 
}
\usage{
tacsatMCP(tacsat, pctThreshold=90)
}
\arguments{
  \item{tacsat}{tacsat data frame}
  \item{pctThreshold}{Percentage of pings to be included into the MCP}
}
\details{
This function is based on the calc_mcp function of the aspace package. The plot_mcp function allows to plot the polygon.
It is used in the DCF Indicator 6 as a filtering of aggregated fisheries.
}
\value{
The tacsat is returned with attached a column called INMCP which flag the pings inside the MCP.
The nodes of the MCP are also returned.
}
\references{EU lot 2 project}
\author{Fabrizio Manco}
\seealso{\code{\link{indicators}}}
\examples{
require(PBSmapping)
data(tacsat)

# Flag the pings inside a polygon gathering 90\% of the pings
tacsat<-tacsatMCP(tacsat, pctThreshold=90)
# Filter the tacsat data to remove the points outside the MCP
tacsat<-subset(tacsat, tacsat$INMCP!=0)
# Grid the filtered vms points
vmsgrid<-vmsGridCreate(tacsat, nameLon = "SI_LONG", nameLat = "SI_LATI",
                       cellsizeX=0.05, cellsizeY=0.05, plotMap=TRUE)
# Add the MCP to the map
plot_mcp(plotnew=FALSE, plotpoints=FALSE, titletxt="")
}
