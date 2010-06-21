\name{findEndVMS}
\alias{findEndVMS}
\title{Finding the succeeding VMS datapoint based on an interval with a specified margin}
\description{
To create an interpolation, two succeeding VMS datapoints are needed. This function
finds the succeeding VMS point and tests if the point is within the specified
time interval and margins. As well, if no succeeding datapoint can be found, this
information is returned
}
\usage{
findEndVMS(VMS, startVMS, interval, margin)
Also used within 'interpolateVMS()'.
}
\arguments{
  \item{VMS}{The VMS dataset}
  \item{startVMS}{Index of VMS dataset of startpoint of interpolation}
  \item{interval}{Time in minutes between the succeeding datapoints}
  \item{margin}{Deviation from specified interval to find succeeding datapoints}
}
\details{
Interval: In most VMS datasets the succeeding datapoint can be found 1 or 2 hours appart. This interval time
should be specified here. Interval can also be specified as e.g. 15 minutes if the VMS / GPS dataset allows this.
Margin: Hardly ever, the interval time is precise. To allow some deviation from the 1 or 2 hour interval
the margin can be adjusted.

The result returned consists of 2 values. The first value is the index of the VMS set specified of the
succeeding datapoint. The second value indicates if the dataset has ended.
If 1st: NA and 2nd 0 then no succeeding VMS point could be found in the specified interval.
If 1st: NA and 2nd 1 then no succeeding VMS point could be found and end of dataset for a specific vessel has been reached.
If 1st: NA and 2nd 2 then no succeeding VMS point could be found and end of complete dataset has been reached.
If 1st: value then 2nd will be 0, succeeding VMS point is found and is specified in 1st value.
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\note{
This function is called inside interpolateVMS()
}
\seealso{\code{filterVMS(),interpolateVMS()}}
\examples{
data(VMS)
startVMS <- 2
findEndVMS(VMS,startVMS,interval=120,margin=10) #result: 3 0 Succeeding point = VMS[3,] and end dataset has not been reached yet.
}

