\name{segmentedTacsatSpeed}
\alias{segmentedTacsatSpeed}
\title{
Define activity based on segmented regression of speed profile}
\description{
Given the speed profile by gear or vessel in a user defined time frame a segmented
regression analyses will be performed to indicate where the fishing speeds are
located. The segmented regression takes place on the cumulative speed profile.
}
\usage{
segmentedTacsatSpeed(tacsat, units = "year", analyse.by = "VE_REF",
                     speed = "calculated", logfit = FALSE, CI = 0.95,
                     saveDir = tempdir())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tacsat}{
A tacsat dataset (with optional column "LE_GEAR" when matched to eflalo)
}
  \item{units}{
Analyse by: "year", "month" and "week". "month" and "week"
cannot be used at same time.
}
  \item{analyse.by}{
Analyse tacsat by gear ("LE_GEAR"), vessel ("VE_REF") or a combination of gear
and vessel ("VE_REF+LE_GEAR").
}
  \item{speed}{
Define if speed profile used must be taken as given in the tacsat file
 (speed = "instantanious") or if internally speed must be calculated
(speed = "calculated"). Default is "calculated"
}
  \item{logfit}{
Logical. Define whether the speed profile frequencies must be log-transformed.
Default is F.
}
  \item{CI}{
Define confidence interval for calculated segmented break points. Default is
0.95.
}
  \item{saveDir}{
Directory to save overview and success of fit. Default = tempdir().
}
  \item{forceLowerBound}{
  Fix the lower breakpoint value at the forceLowerBound given. If not specified,
  lower breakpoint is estimated.
}
}
\details{
To fit a speed profile at least 20 VMS pings must exist. The function highly
depends on accurate starting points. After 20 random tries the fitting procedure
exits and returns a '0' success.
}
\value{
SI_STATE = nf for no-fishing and SI_STATE = f for fishing}
\references{
Bastardie et al. 2010
}
\author{
Niels T. Hintzen, Francois Bastardie
}
\seealso{
\code{\link{activityTacsatAnalyse}}, \code{\link{activityTacsat}}

}
\examples{
data(tacsat)
tacsat <- tacsat[1:20000,]

#-Fit based on vessel and calculated speed
newTacsat <- segmentedTacsatSpeed(tacsat,units="year",analyse.by="VE_REF",
                                  speed="calculated",logfit=FALSE,CI=0.95)

data(eflalo)
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]

#-Fit based on gear and instantaneous speed
newTacsat <- segmentedTacsatSpeed(tacsatp,units="year",analyse.by="LE_GEAR",
                                  speed="instantaneous",logfit=FALSE,CI=0.95)

}