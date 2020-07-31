\name{VMShf}
\alias{VMShf}
\docType{data}
\title{Small high frequency test VMS dataset CANNOT BE DISTRIBUTED WITHOUT PRIOR PERMISSION -> SEE Author}
\description{
A short VMS dataset with high frequency interval rate (every 6 minutes) used to test spatial and VMS functions.
Dataframe contains vessel ID, decimal longitude and decimal latitude position, speed and heading at ping
and date of ping.
}
\usage{data(VMShf)}
\format{
  A data frame with 702 observations on the following 6 variables.
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
niels.hintzen@wur.nl
}
\references{
niels.hintzen@wur.nl
}
\examples{
data(VMShf)
str(VMShf)

}
