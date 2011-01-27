\name{clipPolygons}
\alias{clipPolygons}
\title{
Function to clip Polygons together
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clipPolygons(shapeAll, europa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shapeAll}{
%%     ~~Describe \code{shapeAll} here~~
}
  \item{europa}{
%%     ~~Describe \code{europa} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (shapeAll, europa) 
{
    if (!exists("europa")) {
        load("europa.rda")
    }
    allSourcePoly <- numeric()
    for (x in 1:length(shapeAll$shp$shp)) {
        eachPolyFromShape <- shapeAll$shp$shp[[x]]
        sourcePoly <- cbind(eachPolyFromShape$record, seq(1:length(eachPolyFromShape$points$X)), 
            eachPolyFromShape$points)
        allSourcePoly <- rbind(allSourcePoly, sourcePoly)
    }
    colnames(allSourcePoly) <- c("PID", "POS", "X", "Y")
    allPolyFromShape <- as.PolySet(allSourcePoly)
    clipShapeFromLand <- joinPolys(allPolyFromShape, europa, 
        operation = "DIFF", maxVert = 1e+06)
    clipShapeFromLand <- combinePolys(clipShapeFromLand)
    clipShapeFromLand <- cbind(clipShapeFromLand$SID, clipShapeFromLand$POS, 
        clipShapeFromLand$X, clipShapeFromLand$Y)
    colnames(clipShapeFromLand) <- c("PID", "POS", "X", "Y")
    clipShapeFromLand <- as.PolySet(clipShapeFromLand)
    return(clipShapeFromLand)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
