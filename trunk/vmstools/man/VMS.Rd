\name{VMS}
\alias{VMS}
\docType{data}
\title{Small test VMS dataset CANNOT BE DISTRIBUTED WITHOUT PRIOR PERMISSION -> SEE Author}
\description{
A short VMS dataset used to test spatial and VMS functions.
Dataframe contains vessel ID, decimal longitude and decimal latitude position, speed and heading at ping
and date of ping.
}
\usage{data(VMS)}
\format{
  A data frame with 58190 observations on the following 6 variables.
  \describe{
    \item{\code{ship}}{ID of ship}
    \item{\code{declat}}{Decimal notation of longitude position}
    \item{\code{declon}}{Decimal notation of latitude postion}
    \item{\code{speed}}{Speed at ping}
    \item{\code{heading}}{Heading at ping}
    \item{\code{date}}{Date and time of ping}
  }
}
\source{
niels.hintzen 'at' wur.nl
}
\references{
niels.hintzen 'at' wur.nl
}
\examples{
data(VMS)
str(VMS)

}

