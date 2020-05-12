\name{calc_mcp}
\alias{calc_mcp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing the Minimum Convex Polygon (MCP)
}
\description{
This function computes the Minimum Convex Polygon (MCP) from a set of points. The MCP is the minimum area polygon containing a set of point locations.
}
\usage{
calc_mcp(id=1, points = activities, filename="MCP_Output.txt",
verbose = FALSE, pct = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
Provide a unique integer to identify an MCP from others that you may construct with other data points
}
  \item{points}{
Two-column matrix or data frame of point coordinates
}
  \item{filename}{
A character name for an ASCII output file
}
  \item{verbose}{
Boolean: set to TRUE if extended processing feedback is wanted
}
  \item{pct}{
Integer 0 <= pct <=100, the percentage of the MCP for which area is provided
}
}
\details{
This function is most powerful when used repetitively within a loop to compute the MCP for subsets of points stored in a large data table.
}
\value{
The returned result is a list:
%%  If it is a LIST, use
\item{MCP.area }{The area of the MCP in square kilometers}
\item{MCP.pct }{The desired percentage of the MCP for which the area is computed}
\item{MPC.coords}{A matrix containing MCP vertices. Each row represents a unique point,
the first column contains x-coordinates, and the second, y-coordinates
}}
\references{
Builds upon MCP functions available in the adehabitat package
}
\author{
Randy Bui, Ron N. Buliung, Tarmo K. Remmel
}
\note{
Results are stored in the r.MCP object (required for graphical visualization using
plot_mcp). This function can be used on its own (once) or repetitively in a loop
to process grouped point data stored in a larger table. When used repetitively, be
sure to increment the id parameter to ensure that each MCP has a unique identifier.
The output ASCII coordinate file can be further processed using the makeshapes
function to generate an ESRI Shapefile for MCP polygons.}

\examples{
data(tacsat)
calc_mcp(id=1, points = tacsat[1:10,c("SI_LONG","SI_LATI")], filename="MCP_Output.txt",
verbose = FALSE, pct = 100)
}