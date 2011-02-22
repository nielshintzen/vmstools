\name{eflalo}
\alias{eflalo}
\docType{data}
\title{
Vessel logbook dataset
}
\description{
A dataset consisting of logbook (landings / values) registrations of disguised origin for 2 consecutice years.
}
\usage{data(eflalo)}
\format{
  A data frame with 4539 observations on the following 189 variables.
  \describe{
    \item{\code{VE_REF}}{Vessel reference / ID}
    \item{\code{VE_FLT}}{Fleet reference (DCF regulation)}
    \item{\code{VE_COU}}{Flag nation of vessel}
    \item{\code{VE_LEN}}{Vessel length}
    \item{\code{VE_KW}}{Vessel power (in kW)}
    \item{\code{VE_TON}}{Vessel Tonnage (GT)}
    \item{\code{FT_REF}}{Fishing trip reference number}
    \item{\code{FT_DCOU}}{Departure country}
    \item{\code{FT_DHAR}}{Departure harbour}
    \item{\code{FT_DDAT}}{Departure date}
    \item{\code{FT_DTIME}}{Departure time}
    \item{\code{FT_LCOU}}{Landing country}
    \item{\code{FT_LHAR}}{Landing harbour}
    \item{\code{FT_LDAT}}{Landing date}
    \item{\code{FT_LTIME}}{Landing time}
    \item{\code{LE_ID}}{Log event ID}
    \item{\code{LE_CDAT}}{Catch date}
    \item{\code{LE_STIME}}{Log event start time}
    \item{\code{LE_ETIME}}{Log event end time}
    \item{\code{LE_SLAT}}{Log event start position latitude}
    \item{\code{LE_SLON}}{Log event start position longitude}
    \item{\code{LE_ELAT}}{Log event end position latitude}
    \item{\code{LE_ELON}}{Log event end position longitude}
    \item{\code{LE_GEAR}}{Gear}
    \item{\code{LE_MSZ}}{Mesh size}
    \item{\code{LE_RECT}}{ICES rectangle}
    \item{\code{LE_DIV}}{ICES division}
    \item{\code{LE_MET_level6}}{Fishing activity}
    \item{\code{LE_KG_SPEC}}{Landing weight estimate of species (FAO species codes)}
    \item{\code{LE_EURO_SPEC}}{Landing value of species (FAO species codes)}
   }
}
\details{
This example dataset if for example purposes only. Without prior permission of the authors it is not allowed to use this data
other than for example non-publishable purposes.
}
\source{
NA
}
\references{
EU Lot 2 project, Niels T. Hintzen (niels.hintzen@wur.nl)
}
\examples{
data(eflalo)
eflalo <- formatEflalo(eflalo)
str(eflalo)
}
\keyword{datasets}
