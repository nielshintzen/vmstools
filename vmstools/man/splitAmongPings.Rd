\name{splitAmongPings}
\alias{splitAmongPings}
\title{
Split values or landings from eflalo over tacsat pings
}
\description{
Split the values or landings as listed in the eflalo file over the tacsat pings,
while taking different levels into account such as by day, ICESrectangle or by trip number.
Also there is a possibility to merge the eflalo records without a matching tacsat trip.
}
\usage{
splitAmongPings(tacsat, eflalo, variable = "all",
                level = "day", conserve = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tacsat}{
Tacsat object
}
  \item{eflalo}{
Eflalo object
}
  \item{variable}{
Indicating what to split: "all","value","kgs"
}
  \item{level}{
Levels can be: "day", "ICESrectangle", "trip"
}
  \item{conserve}{
Logical, if kgs or value needs to be conserved if merging by trip number is not possible (default = TRUE)
}
}
\details{
Levels have hierachy, so if "day" is specified, also "ICESrectangle" and "trip" will be used.
If "ICESrectangle" is specified also "trip" will be used.
"Trip" can be used on its own. Same hierachy applies to merging when conserve = TRUE (except for trip level).

Note that tacsat file needs a column SI_STATE which has 0 for non-fishing and 1 for fishing records.
}
\value{
Merged tacsat file will be returned including the splitted values over the tacsat pings where SI_STATE is not zero.
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{
\code{\link{mergeEflalo2Tacsat}}, \code{\link{mergeEflalo2Pings}}}
\examples{
data(tacsat); tacsat <- tacsat[1:1000,]
data(eflalo); eflalo <- eflalo[1:1000,]

tacsatp           <- mergeEflalo2Tacsat(eflalo,tacsat)

#- Create a column names SI_STATE which holds values 0 or 1 which denotes no
#  fishing & fishing.
tacsatp$IDX       <- 1:nrow(tacsatp)
tacsatFilter      <- filterTacsat(tacsatp,st=c(1,6),hd=NULL,remDup=TRUE)
tacsatp$SI_STATE  <- 0
tacsatp$SI_STATE[tacsatFilter$IDX] <- 1

tacsatp           <- subset(tacsatp,SI_STATE == 1)

tacsatEflalo      <- splitAmongPings(tacsat=tacsatp,eflalo=eflalo,
                        variable="all",level="day",conserve=TRUE)
}