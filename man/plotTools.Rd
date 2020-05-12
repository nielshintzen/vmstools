\name{plotTools}
\alias{plotTools}
\title{
Plot eflalo or tacsat files
}
\description{
Plot eflalo or tacsat given a grid and column name as a map of intensity
}
\usage{
plotTools(x,level= "ICESrectangle",xlim,ylim,zlim=NULL,log=FALSE,
          gridcell=c(0.1,0.05),color=NULL,control.tacsat= list(clm = NULL),
          control.eflalo = list(clm = NULL), returnRange = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Eflalo or tacsat dataframe
}
  \item{level}{
  Aggregating level: "ICESrectangle" or specified "gridcell" between xlim and ylim in
  steps gricell for tacsat or "ICESrectangle" only for eflalo
}
  \item{xlim}{
  two element numeric vector giving a range of longitudes, expressed in degrees,
  to which drawing should be restricted. Longitude is measured in degrees east
  of Greenwich, so that, in particular, locations in the USA have negative
  longitude. If fill = TRUE, polygons selected by region must be entirely inside
   the xlim range. The default value of this argument spans the entire longitude
   range of the database.
}
  \item{ylim}{
  two element numeric vector giving a range of latitudes, expressed in degrees,
  to which drawing should be restricted. Latitude is measured in degrees north
  of the equator, so that, in particular, locations in the USA have positive
  latitude. If fill = TRUE, polygons selected by region must be entirely inside
  the ylim range. The default value of this argument spans the entire latitude
  range of the database.
}
  \item{zlim}{
  the minimum and maximum z values for which colors should be plotted,
  defaulting to the range of the finite values of z. Each of the given colors
  will be used to color an equispaced interval of this range. The midpoints of
  the intervals cover the range, so that values just outside the range will be
  plotted.
}
  \item{log}{
  whether values to plot need to be logged, TRUE or FALSE
}
  \item{gridcell}{
  two element numeric vector giving the steps in longitudinal degrees and steps
  in latitudinal degrees for the grid to plot on. Only needed when
  level = "gridcell"
}
  \item{color}{
  colors range. default = brewer.pal(9,"YlOrRd")
}
  \item{control.tacsat}{
  list with attribute "clm": column names to display aggregated results over
}
  \item{control.eflalo}{
  list with attribute "clm": column names to display aggregated results over
}
  \item{returnRange}{Logical: return range of plotted values (which can be used to
  define own legend}
  
  \item{las}{Direction of axis legends. See ?par for more info}
}
\author{Niels T. Hintzen}
\seealso{
\code{\link{plotTreeMap}}, \code{\link{Grid2KLM}}, \code{\link{landingsMap2GIFanim}}, \code{\link{pings2EffortMaps}}, \code{\link{pings2LandingsMaps}}}

\examples{
data(tacsat)
data(eflalo)

plotTools(tacsat,level="ICESrectangle",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,
          log=FALSE,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm=NULL))
x11()
plotTools(eflalo,level="ICESrectangle",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,
          log=FALSE,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm=NULL))
x11()
plotTools(tacsat,level="gridcell",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,
          log=FALSE,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm=NULL))

x11()
plotTools(eflalo,level="ICESrectangle",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,
          log=FALSE,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm=NULL),
          control.eflalo=list(clm=c("LE_KG_COD","LE_KG_PLE")))

}
