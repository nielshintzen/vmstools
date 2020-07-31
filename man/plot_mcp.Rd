\name{plot_mcp}
\alias{plot_mcp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Minimum Convex Polygon
}
\description{
This function plots the MCP as a polygon, which covers the geographical extent of a set of points on a Cartesian plane.
}
\usage{
plot_mcp(plotnew=TRUE, plotpoints=TRUE, points.col='black',
points.pch=1, titletxt="Title", xaxis= "Easting (m)",
yaxis="Northing (m)", mcp.col='black', mcp.lwd=2,
fill.col=NA, jpeg=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plotnew}{
Boolean: Set to TRUE to create a new plot. Set to FALSE to overlay current plot.}
  \item{plotpoints}{
Boolean: Set to TRUE if the point observations are to be plotted}
  \item{points.col}{
Specify a colour for the point observations}
  \item{points.pch}{
Specify a plotting symbol for the point observations}
  \item{titletxt}{
A string to use as the title on the plot}
  \item{xaxis}{
A string to label the x-axis of the plot
}
  \item{yaxis}{
A string to label the y-axis of the plot
}
  \item{mcp.col}{
Specify the line colour for the MCP
}
  \item{mcp.lwd}{
Specify the line width for the MCP
}
  \item{fill.col}{
Specify a fill colour for the MCP
}
  \item{jpeg}{
Boolean: Set to TRUE if the plot should be saved in JPEG format
}
  \item{\dots}{
Arguments to be passed to graphical parameters
}
}
\details{
The r.MCP object (generated in calc_mcp function) is required to plot the MCP.
}
\author{
Randy Bui, Ron N. Buliung, Tarmo K. Remmel
}

\examples{
data(tacsat)

calc_mcp(id=1, points = tacsat[1:10,c("SI_LONG","SI_LATI")], filename="MCP_Output.txt",
verbose = FALSE, pct = 100)
plot_mcp(plotnew=TRUE, plotpoints=TRUE, titletxt="Title",
xaxis= "Easting (m)", yaxis="Northing (m)")
}