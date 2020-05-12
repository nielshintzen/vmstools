\name{countPings}
\alias{countPings}
\title{
Count the number of VMS pings in a selection
}
\description{
Count the number of VMS pings in any selection made in time and spatial frame.
}
\usage{
countPings(formula, tacsat, grid = NULL)
}
\arguments{
  \item{formula}{ specify the elements you want to use as axis of the count}
  \item{tacsat}{ tacsat data file, possibly with additional columns}
  \item{grid}{ if in formulate 'gridcell' is chosen, a SpatialGrid must be provided}
}
\details{
Formula has form: ~Variable+timeVariable+anotherTimeVariable+spatialVariable \cr

options in formula for Variable: any column name of tacsat file   \cr
options in formula for timeVariable: day, week, month, quarter, year    \cr
options in formula for spatialVariable: icesrectangle, icesarea, gridcell \cr
}
\value{
Returns the matrix with counted pings by each specified variable
}
\references{Hintzen et al. 2010 Improved estimation of trawling tracks using cubic Hermite spline interpolation of position registration data,
EU lot 2 project}
\author{Niels T. Hintzen}
\note{
if Tacsat is a big file, the overlay function might fail resulting in terminating the function}
\seealso{
\code{\link{createGrid}}, \code{\link{vmsGridCreate}}
}
\examples{
data(tacsat)

#make the tacsat file a bit smaller
tacsat  <- tacsat[1:10000,]

grid          <- createGrid(range(tacsat$SI_LONG,na.rm=TRUE),
                  range(tacsat$SI_LATI,na.rm=TRUE),0.5,0.5,type="SpatialGrid")

result        <- countPings(~VE_REF+year+gridcell,tacsat,grid=grid)
result        <- countPings(~VE_REF+week+year+icesrectangle,tacsat)
}

