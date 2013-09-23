\name{clipObs2Tacsat}
\alias{clipObs2Tacsat}
\title{
Find tacsat records close to observation data records
}
\description{
Find the tacsat records that are close in distance or on a defined grid to
the observation data positions in both time and space
}
\usage{
clipObs2Tacsat(tacsat, obs, method, control.grid, control.euclidean,
               temporalRange, all.t = FALSE, rowSize)
}
\arguments{
  \item{tacsat}{
    Dataframe of vessel VMS in tacsat format
}
  \item{obs}{
    Dataframe also in tacsat format, but with observation data (GPS and Date+Time)
}
  \item{method}{
  'Grid' or 'Euclidian'. Match tacsat records in same gridcell (method = "grid")
  or only x km from observed datapoint (method = "euclidean")
}
  \item{control.grid}{
   If method = "grid", control.grid is a list with specification of the spatialGrid
   (spatGrid, of class "SpatialGrid") or the resolution in x and y values where
   the outermargins can be given by the tacsat dataset (gridBox="tacsat") or the
   observation dataset (gridBox="obs").

   Example: control.grid=list(spatGrid = NULL, resx = c(0.5,1),
                              resy = NULL, gridBbox = "obs")
}
  \item{control.euclidean}{
   If method = "euclidean" then a maximum distance the tacsat position may deviate
   from the observed datapoint needs to be given (in km).
}
  \item{temporalRange}{
  The range in which tacsat records may deviate from the observation time stamp (in minutes)
}
  \item{all.t}{
  If you want to return all tacsat records, all.t=T (return = list, first
  element = observation set, second = matching tacsat), else F second = all tacsat
}
  \item{rowSize}{
  To speed up calculations, define maximum rowSize (default = 1000)
}
}
\details{
Returned is a list with the first element being the observation dataset with an
GR_ID identifier that matches the GR_ID from the matching tacsat records. The
second element is the tacsat records with the same GR_IDs.
}
\references{EU Lot 2 project
}
\author{
Niels T. Hintzen
}
\examples{
data(tacsat)
tacsat        <- sortTacsat(tacsat)

obs           <- tacsat[round(seq(1,nrow(tacsat),length.out=30)),]
obs           <- obs[,c(1,2,3,4,5,6)]
obs$OB_TYP    <- "Benthos"
colnames(obs) <- c("OB_COU","OB_REF","SI_LATI","SI_LONG","SI_DATE","SI_TIME","OB_TYP")
newTime       <- obs$SI_DATIM - runif(30,-60*20,60*20)
obs$SI_LATI   <- jitter(obs$SI_LATI,factor=0.25)
obs$SI_LONG   <- jitter(obs$SI_LONG,factor=0.5)


a <- clipObs2Tacsat(tacsat,obs,method="grid",control.grid=list(resx=0.1,resy=0.05,
                    gridBbox="obs"),temporalRange=c(-30,120),all.t=FALSE)
a <- clipObs2Tacsat(tacsat,obs,method="euclidean",control.euclidean=list(threshold=1),
                    temporalRange=c(-1e10,-1) ,all.t=FALSE)
}

