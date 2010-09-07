\name{findEndTacsat}
\alias{findEndTacsat}
\title{Finding the succeeding Tacsat datapoint based on an interval with a specified margin}
\description{
To create an interpolation, two succeeding Tacsat datapoints are needed. This function
finds the succeeding Tacsat point and tests if the point is within the specified
time interval and margins. As well, if no succeeding datapoint can be found, this
information is returned
}
\usage{
findEndTacsat(tacsat, startTacsat, interval, margin)
Also used within 'interpolateVMS()'.
}
\arguments{
  \item{tacsat}{The Tacsat dataset}
  \item{startTacsat}{Index of Tacsat dataset of startpoint of interpolation}
  \item{interval}{Time in minutes between the succeeding datapoints}
  \item{margin}{Deviation from specified interval to find succeeding datapoints}
}
\details{
Interval: In most Tacsat datasets the succeeding datapoint can be found 1 or 2 hours appart. This interval time
should be specified here. Interval can also be specified as e.g. 15 minutes if the Tacsat / GPS dataset allows this.
Margin: Hardly ever, the interval time is precise. To allow some deviation from the 1 or 2 hour interval
the margin can be adjusted.

The result returned consists of 2 values. The first value is the index of the Tacsat set specified of the
succeeding datapoint. The second value indicates if the dataset has ended.
If 1st: NA and 2nd 0 then no succeeding Tacsat point could be found in the specified interval.
If 1st: NA and 2nd 1 then no succeeding Tacsat point could be found and end of dataset for a specific vessel has been reached.
If 1st: NA and 2nd 2 then no succeeding Tacsat point could be found and end of complete dataset has been reached.
If 1st: value then 2nd will be 0, succeeding Tacsat point is found and is specified in 1st value.
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\note{
This function is called inside interpolateTacsat()
}
\seealso{\code{filterTacsat(),interpolateTacsat()}}
\examples{
data(tacsat)
startTacsat <- 2
findEndTacsat(tacsat,startTacsat,interval=120,margin=10) #result: 3 0 Succeeding point = tacsat[3,] and end dataset has not been reached yet.
}

