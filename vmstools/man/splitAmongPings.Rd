\name{splitAmongPings}
\alias{splitAmongPings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splitAmongPings(tacsat, eflalo, variable = "all", level = "day", conserve = T)
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
Logical, if kgs or value needs to be conserved if merging by trip number is not possible (default = T)
}
}
\details{
Levels have hierachy, so if "day" is specified, also "ICESrectangle" and "trip" will be used.
If "ICESrectangle" is specified also "trip" will be used.
"Trip" can be used on its own. Same hierachy applies to merging when conserve = T (except for trip level).
}
\value{
Merged tacsat file will be returned including the splitted values over the tacsat pings where SI_STATE = fishing.
}
\references{EU Lot 2 project}
\author{Niels T. Hintzen}
\seealso{
\code{mergeEflalo2Tacsat(), mergeEflalo2Pings()}
}
\examples{
data(tacsat)
data(eflalo)

tacsatp           <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatEflalo      <- splitAmongPings(tacsat=tacsatp,eflalo=eflalo,variable="all",level="day",conserve=T)
}