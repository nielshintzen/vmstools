\name{effort}
\alias{effort}
\title{
Calculate effort of tacsat or eflalo dataset
}
\description{
Calculate effort (in hours, days, ...) based on tacsat or eflalo dataset for different combination of settings}
\usage{
effort(x, by = "FT_REF", unit = "hours", weight = c(0.5, 0.5), fill.na = FALSE)
}
\arguments{
  \item{x}{
Either eflalo or tacsat format data}
  \item{by}{
Vector including the elements to calculate effort by. Options are: byRow, VE_REF, FT_REF, LE_GEAR, SI_DAY, SI_WEEK, SI_MONTH, SI_QUARTER, SI_YEAR, LE_RECT, LE_AREA
}
  \item{unit}{
Unit must be in 'secs,mins,hours,days or weeks'}
  \item{weight}{
Only relevant for tacsat: weight to apply to calculation of mean interval rate towards and away from ping
}
  \item{fill.na}{
Only relevant for tacsat: If interval rate cannot be calculated based on default or provided weight, take closest alternative to provide an interval rate
}
}
\details{
if 'byRow' is selected, no other elements can be added}
\author{
Niels T. Hintzen}
\seealso{
\code{\link{raiseTacsat}}
}
\examples{
data(eflalo)
data(tacsat)
#-Remove duplicated records from tacsat
myf       <- paste(tacsat$VE_REF,tacsat$SI_LATI,tacsat$SI_LONG,
                   tacsat$SI_DATE,tacsat$SI_TIME);
tacsat    <- tacsat[!duplicated(myf),];

#- Try some out for eflalo
a1 <- effort(eflalo,by=c("FT_REF"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(a1$EFFORT)
a2 <- effort(eflalo,by=c("SI_MONTH"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(a2$EFFORT)
a3 <- effort(eflalo,by=c("VE_REF","SI_MONTH"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(a3$EFFORT)
a4 <- effort(eflalo,by=c("VE_REF","SI_QUARTER","LE_GEAR"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(a4$EFFORT)
a5 <- effort(eflalo,by=c("byRow"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(a5$EFFORT)
a6 <- effort(eflalo,by=c("SI_DAY","LE_GEAR"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(a6$EFFORT)

#- Try some out for tacsat
tacsatp           <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp$LE_GEAR   <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]
b1 <- effort(tacsatp,by=c("FT_REF"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(b1$EFFORT)
b2 <- effort(tacsatp,by=c("SI_MONTH"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(b2$EFFORT)
b3 <- effort(tacsatp,by=c("VE_REF","SI_MONTH"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(b3$EFFORT)
b4 <- effort(tacsatp,by=c("VE_REF","SI_QUARTER","LE_GEAR"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(b4$EFFORT)
b5 <- effort(tacsatp,by=c("byRow"),unit="hours",weight=c(0.5,0.5),fill.na=TRUE); sum(b5$EFFORT,na.rm=T)}