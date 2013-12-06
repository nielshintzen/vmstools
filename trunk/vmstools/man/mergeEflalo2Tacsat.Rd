\name{mergeEflalo2Tacsat}
\alias{mergeEflalo2Tacsat}
\docType{data}
\title{
Merge TACSAT formatted data to EFLALO2 formatted data
}
\description{
Method to assign the EFLALO2 tripnumber to the accompanying TACSAT records, which were registered between the departure date and arrival date of a trip.
}
\usage{mergeEflalo2Tacsat(eflalo2,tacsat)}

\arguments{
  \item{eflalo2}{Trip and landings data in the Eflalo2 format}
  \item{tacsat}{Vms data in the Tacsat format}
}
\details{
Method only assignes a tripnumber to the tacsat data if there are vms pings that occur between the departure date and arrival date of a specific vessel. If no 
related trip can be found, tripnumber assigned to tacsat equals NA.
}
\references{EU lot 2 project }
\author{Niels T. Hintzen, Francois Bastardie}
\seealso{\code{\link{findEndTacsat}}, \code{\link{sortTacsat}},
               \code{\link{filterTacsat}}}
\examples{
  #Load the eflalo2 and tacsat data
data(eflalo)
data(tacsat)

#-Remove duplicated records from tacsat
myf       <- paste(tacsat$VE_REF,tacsat$SI_LATI,tacsat$SI_LONG,
                   tacsat$SI_DATE,tacsat$SI_TIME);
tacsat    <- tacsat[!duplicated(myf),];

#-Merge the tacsat and eflalo2 data
#computation time approx 30 seconds
result <- mergeEflalo2Tacsat(eflalo,tacsat)
}
