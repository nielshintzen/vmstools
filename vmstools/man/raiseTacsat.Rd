\name{raiseTacsat}
\alias{raiseTacsat}
\title{
Raise the effort in tacsat to the effort available in eflalo.
}
\description{
The eflalo dataset contains a total amount of effort. This total amount can be
set as the total amount of effort to be found in tacsat too. This function raises
the total amount of effort in tacsat to that of eflalo.
}
\usage{
raiseTacsat(tacsat, eflalo, by = c("LE_GEAR", "VE_REF", "SI_DAY", "LE_RECT", "SI_YEAR"), sortBy = TRUE)
}
\arguments{
  \item{tacsat}{
Tacsat object
}
  \item{eflalo}{
Eflalo object
}
  \item{by}{
By can be: "SI_DAY", "SI_WEEK", "SI_MONTH", "SI_YEAR","LE_RECT","LE_ICESAREA","SI_QUARTER",
}
  \item{sortBy}{
Logical. If items in 'by' should be ordered by descending uniqueness (if TRUE match between tacsat and eflalo increases)}
}
\details{
The function looks also for a 'PropFish' column. This allows users to scale down
the total effort in eflalo according to the proportion that is assumed to be fishing
effort. This helps to raise the tacsat data where only records which have SI_STATE = 1
will be used}

\author{
Niels T. Hintzen}
\seealso{
\code{\link{splitAmongPings}}, \code{\link{estimatePropFishing}}}
\examples{
data(tacsat)
data(eflalo)

#- Calculate effort in tacsat dataset and define activity (fishing)
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp <- sortTacsat(tacsatp)

#- Define effort on a trip level
tacsatp <- intervalTacsat(tacsatp,level="trip",fill.na=T,weight=c(0.5,0.5))
tacsatp1<- subset(tacsatp,FT_REF==0)
tacsatp2<- subset(tacsatp,FT_REF!=0)
#- Define effort on a vessel level
tacsatp1<- intervalTacsat(tacsatp1,level="vessel",fill.na=T,weight=c(0.5,0.5))

#- Combine the two datasets with effort associated
tacsatp <- rbind(tacsatp1,tacsatp2)
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]

#- Quick and dirty approach to defin activity
tacsatp <- filterTacsat(tacsatp,st=c(4,8))
tacsatp$INTV[which(tacsatp$INTV>232)] <- 115
tacsatp$INTV[which(is.na(tacsatp$INTV)==T)] <- 115
tacsat  <- tacsatp
tacsat$SI_STATE <- 1

#- Estimate the proportion fishing first (not necessary)
eflalo  <- estimatePropFishing(tacsat,eflalo,by=c("VE_REF","FT_REF"))

#- Raise the tacsat dataset to the effort of the eflalo dataset
tacsat  <- raiseTacsat(tacsat,eflalo,by=c("LE_GEAR","SI_DAY","SI_YEAR"),sortBy=T)
tacsat  <- raiseTacsat(tacsat,eflalo,by=c("LE_GEAR","SI_DAY","SI_YEAR"),sortBy=F)

}
