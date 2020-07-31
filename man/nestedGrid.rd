\name{nestedGrid}
\alias{nestedGrid}
\title{
Define a nested grid}
\description{
Assign VMS points to a nested grid so that areas with a high density of points
will have small grid cells and areas with a low density will have larger cells.
}
\usage{
nestedGrid(tacsat, resx, resy, maxstep = 10, n = 20, control = list(clm = NULL, FUN = NULL))
}
\arguments{
  \item{tacsat}{
Tacsat dataframe
}
  \item{resx}{
gridcell size in degrees in the longitude / x direction
}
  \item{resy}{
gridcell size in degrees in the latitude / y direction
}
  \item{maxstep}{
The maxiumum number of times the grid cells should be split
}
  \item{n}{
If the number of points in a cell \code{>= n} then split the cell
}
  \item{control}{
A list determining the output in the \code{data} slot, the possible components are:
\describe{
  \item{\code{clm}}{
    The name of the tacsat column that \code{FUN} should be applied to.
    Also the name of the only column in the dataframe that makes up the \code{data} slot.
    The default (\code{clm = NULL}) results \code{FUN} being applied a new
    column called \code{"count"} which is simply \code{rep(1,nrow(tacsat))}.
  }
  \item{\code{FUN}}{
    The function to be applied to \code{tacsat[[clm]]}.
    The default (\code{FUN = NULL}) results in the function \code{sum}, so
    if \code{clm = NULL} and \code{FUN = NULL} the result will be a count
    of the number of datapoints.
  }
}
}

}
\details{
The alogrithm works as follows: A coarse starting grid is applied to the positional
data, the number of datapoints in each grid cell is counted and if this number \code{>= n}
then the cell is split in two.  Now the number of datapoints in the smaller cells are
counted again and any cells with \code{>= n} will be split again, up to \code{maxstep}
times.

This function allows data-rich areas to be plotted with a high spatial resolution while
allowing a lower spatial resolution for data-poor areas.  The nested grid also
tends to reduce the amount of clustering within each grid cell, which is important
for estimating the area impacted by fishing gear inside each cell.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A SpatialPolygonsDataFrame, the value of the \code{data} slot will be
determined by the \code{control} settings.
}
\references{
Gerritsen, H. D., Minto, C. and Lordan, C. (2013) How much of the seabed is
impacted by mobile fishing gear? Absolute estimates from Vessel Monitoring
System (VMS) point data. ICES Journal of Marine Science \bold{XX:XX},
doi: 10.1093/icesjms/fst017
}
\author{
Hans D. Gerritsen, Niels T. Hintzen

}
\seealso{
\code{\link{polyDef}}
}
\examples{
data(tacsat)
tacsat            <- tacsat[sample(nrow(tacsat),2500),] # to speed it up a bit
tacsat            <- intervalTacsat(tacsat,level="vessel",fill.na=TRUE)
tacsat$INTV       <- ifelse(tacsat$INTV > 240, 240, tacsat$INTV)
tacsat$GEAR_WIDTH <- 0.024
tacsat$SWEPT_AREA <- tacsat$INTV / 60 * tacsat$SI_SP * tacsat$GEAR_WIDTH
SPDF              <- nestedGrid(tacsat, resx=1, resy=0.5, maxstep=10, n=20,
                      control=list(clm="SWEPT_AREA",FUN=sum))
SP                <- as(SPDF,"SpatialPolygons")
SP                <- surface(SP)
tempfun           <- function(x){lapply(x@Polygons,function(y){y@area})}
SPDF@data$surface <- unlist(lapply(SP@polygons,tempfun))
SPDF@data$SAratio <- SPDF@data$SWEPT_AREA / SPDF@data$surface
breaks            <- c(seq(0,quantile(SPDF@data$SAratio,0.95),length.out=9),2,5,10)
i                 <- cut(SPDF@data$SAratio,breaks=breaks)
SPDF@data$col     <- grey(seq(1, 0, length=length(breaks)))[i]

plot(NA, xlim = c(1, 5), ylim = c(51.5,55), xlab = 'Longitude',ylab = 'Latitude'
  ,asp=1/lonLatRatio(3,53))
plot(SP,col=SPDF@data$col,add=TRUE,border="lightgrey"); box()
points(tacsat$SI_LONG,tacsat$SI_LAT,cex=0.1,col=4)
}