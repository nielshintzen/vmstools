\name{activityTacsat}
\alias{activityTacsat}
\title{Define activity of fishing vessel based on speed profile analyses.
}
\description{
Analyse tacsat data by gear or vessel and define, based on speed profile and
fitting normal distributions through these speed profiles, fishing and other
activities.
}
\usage{
activityTacsat(tacsat, units = "year", analyse.by = "LE_GEAR", storeScheme = NULL,
               plot = FALSE, level="all",sigma=0.911,fixPeaks=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tacsat}{tacsat dataset}
  \item{units}{Analyse by: "year", "month" and "week". "month" and "week" cannot
               be used at same time.}
  \item{analyse.by}{Analyse tacsat by gear ("LE_GEAR") or vessel ("VE_REF").}
  \item{storeScheme}{If \code{\link{activityTacsatAnalyse}} is used, supply
        output here.}
  \item{plot}{Logical. Whether the results of the fit of the normal distributions should be
              plotted}
  \item{level}{If analyse.by="LE_GEAR" is selected, there is an option to analyse
               at vessel level too, but taking the gear parameter estimates as
               starting values. level = "vessel" turns this option on.}
}
\details{
The analysis assumes that the speed profile close to 0 can be mirrored to create
a good normal distribution.
Success of fit might depend on the initial guess in peaks or mean peak values.

Note that the default value of sigma0 = 0.911 (used when prior estimates of parameters
are not given). The value of 0.911 corresponds to a width of 1.5knots on each side
of the mean under a 90 percent CI.
}
\value{
In general, a 5-peak analysis results in h=no fishing / in harbour, f=fishing, s=steaming.
A 3-peak analyses results in f=fishing (closest to zero speed) and s=steaming.
}
\author{Niels T. Hintzen}
\seealso{
\code{\link{activityTacsatAnalyse}}, \code{\link{segmentedTacsatSpeed}}
}
\examples{
data(tacsat)
data(eflalo)

tacsatp         <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]

tacsatp$LE_GEAR <- ac(tacsatp$LE_GEAR)
tacsatp$LE_GEAR[which(is.na(tacsatp$LE_GEAR)==TRUE)] <- "NO_GEAR"
tacsat          <- tacsatp

#--------- LE_GEAR -----------
#- Visual analyses of activity, and define number of peaks / kernels that can be
#  observed (either 3 or 5)
\dontrun{
storeScheme   <- activityTacsatAnalyse(subset(tacsat,LE_GEAR == "OTM" &
                      format(SI_DATIM,"%Y") == 1801),units="year",analyse.by="LE_GEAR",identify="means")
res           <- activityTacsat(subset(tacsat,LE_GEAR == "OTM" &
                      format(SI_DATIM,"%Y") == 1801),units="year",analyse.by="LE_GEAR",
                      storeScheme,plot=TRUE)
res           <- activityTacsat(subset(tacsat,LE_GEAR == "OTM" &
                      format(SI_DATIM,"%Y") == 1801),units="year",analyse.by="LE_GEAR",
                      storeScheme,plot=TRUE,level="vessel")

#--------- VE_REF -----------
tacsat        <- subset(tacsat,VE_REF == "801")
storeScheme   <- activityTacsatAnalyse(tacsat,units="year",analyse.by="VE_REF",
                                       identify="means")

#- Run activityTacsat
res           <- activityTacsat(tacsat,units="year",analyse.by="VE_REF",storeScheme,
                                plot=TRUE,level="all")
}
}
