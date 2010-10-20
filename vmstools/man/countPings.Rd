\name{countPings}
\alias{countPings}
\title{Ping counting over trips}
\description{This function counts the number of pings in the tacsat data for each trip and then merges them with that same file}
\usage{countPings(tacsat = tacsat1f)}
\arguments{
  \item{tacsat}{data.frame of vms data, typically with only fishing pings}
}
\details{}
\references{}
\author{Doug Beare}
\note{}
\seealso{}
\examples{tacsat1f <-countPings(tacsat = tacsat1f) # Where tacsat1f is simply tacsat data with non fishing pings removed
}