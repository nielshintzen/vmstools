\name{Grid2KML}
\alias{Grid2KML}
\title{ Create a grid of colored rectangles for viewing in Google Earth.
}
\description{
The function takes the output from the function vmsGridCreate and scales a png file so it fits
on Google Earth. A WGS84 projection is the default. Data are first converted to a SpatialGridDataFrame
}
\usage{
Grid2KML(output.mat = output.mat, what.quantity = "effort")
}

\arguments{
  \item{output.mat}{This the matrix output by vmsGridCreate}
  \item{what.quantity}{The unit of measurement such as kwhours or hours fished that are displayed in GE and the legend}
}
\details{}
\value{ A png image and a kml file are output. Click on the kml file and the png map will automatically be plotted on GE.

}
\references{ }
\author{Doug Beare
}
\note{}


\seealso{\code{\link{vmsGridCreate}}, \code{\link{createGrid}}}
\examples{

data(tacsat)

grid.size.x <- .05
grid.size.y <- .05

we <- -2
ea <- 9
so <- 48
no <- 58

#Map number of pings

output.mat<-vmsGridCreate(tacsat,we=we,ea=ea,so=so,no=no,
              nameLon = "SI_LONG", nameLat = "SI_LATI",plotMap=TRUE,
              nameVarToSum="",cellsizeX = grid.size.x,
              cellsizeY = grid.size.y,legendx='bottomright',
              numCats=9,plotPoints=F,legendncol=2,legendtitle="npings")


Grid2KML(output.mat=output.mat,what.quantity='npings')

getwd() #Find out where maps are

}
