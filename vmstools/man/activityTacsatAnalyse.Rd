\name{activityTacsatAnalyse}
\alias{activityTacsatAnalyse}
\title{
Analyse speed profiles to define number of peaks.
}
\description{
Analyse histrogram of speed profiles by gear or vessel to define number of peaks
}
\usage{
activityTacsatAnalyse(tacsat, units = "year", analyse.by = "LE_GEAR",
                      identify="peaks")
}
\arguments{
  \item{tacsat}{tacsat dataset}
  \item{units}{Analyse by: "year", "month" and "week". "month" and "week" cannot
               be used at same time}
  \item{analyse.by}{Analyse tacsat by gear ("LE_GEAR") or vessel ("VE_REF")}
  \item{identify}{Identify either the number of peaks ("peaks", rough approach)
                  or mean speed associated with peaks ("means", detailed approach)}
}
\details{
The analyses assumes that the speed profile close to 0 can be mirrored to create
a good normal distribution. The left-hand side of the histogram is identical to the
right-hand side and therefore potential identified peaks should be mirrored as well.
The detailed 'means' approach does take more time to identify but is more succesful
in fitting distributions to the speed profile. It pays off to specify mean values
of the peaks by 0.5knots precision.

Additional Parameters:
Sigma0: Fixed sigma value for distribution around zero. Zero (0) if parameter needs to be estimated.
FixPeaks: Option to fix the top of the density curves at the initial parameter values from 'storeScheme'.

}
\value{
Returns a \code{\link{data.frame}} with a run scheme that can be used as input
to \code{\link{activityTacsat}}
}
\author{Niels T. Hintzen}
\seealso{
\code{\link{activityTacsat}}, \code{\link{segmentedTacsatSpeed}}
}
\examples{
data(tacsat)
data(eflalo)

tacsatp         <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]

tacsatp$LE_GEAR <- ac(tacsatp$LE_GEAR)
tacsatp$LE_GEAR[which(is.na(tacsatp$LE_GEAR)==TRUE)] <- "NO_GEAR"
tacsat          <- tacsatp

\dontrun{
#--------- LE_GEAR -----------
#- Visual analyses of activity, and define number of peaks / kernels that can
#  be observed (either 3 or 5)
storeScheme   <- activityTacsatAnalyse(subset(tacsat,LE_GEAR == "OTM" &
                      year(SI_DATIM) == 1801),units="year",analyse.by="LE_GEAR",
                      identify="peaks")
storeScheme   <- activityTacsatAnalyse(subset(tacsat,LE_GEAR == "OTM" &
                      year(SI_DATIM) == 1801),units="year",analyse.by="LE_GEAR",
                      identify="means")

}
}