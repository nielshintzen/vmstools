\name{nestedGrid}
\alias{nestedGrid}
\title{
Function to create nested grid based on VMS data points}
\description{
Function to create a nested grid based on VMS data points with certain threshold that
decreases the size of each grid cell.
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
Maximum number of grid cell halving steps
}
  \item{n}{
Threshold that when exceeded, grid cell size is halved
}
  \item{control}{
Column to specify over which data needs to be aggregated in function.
FUN specifies which function to use
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
data(tacsat)
tacsat            <- intervalTacsat(tacsat,level="vessel",fill.na=T)
tacsat$GEAR_WIDTH <- 0.024
tacsat$SWEPT_AREA <- tacsat$INTV / 60 * tacsat$SI_SP * tacsat$GEAR_WIDTH

SPDF  <- nestedgrid(tacsat,resx=1,resy=0.5,maxstep=10,n=20,control=list(clm="SWEPT_AREA",FUN=mean))

}