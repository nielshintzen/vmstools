\name{interpolation2Tacsat}
\alias{interpolation2Tacsat}
\title{
Convert the output (a list) from interpolateTacsat to tacsat
}
\description{
This function takes the list of tracks output by interpolateTacsat and converts them back to tacsat format. This means
that all the functions that work on tacsat can work with this new, larger interpolated dataset. For example you could now
distribute the landings among the new pings.  The npoints argument is the optional number of pings between each real or actual ping.
}
\usage{
interpolation2Tacsat(interpolation, tacsat, npoints = 10,equalDist=TRUE)
}

\arguments{
  \item{interpolation}{A list of tracks between actual pings output by the interpolateTacsat function.}
  \item{tacsat}{Tacsat data used when the interpolateTacsat function was called.}
  \item{npoints}{The number of pings or positions required between each real or actual vessel position or ping}
  \item{equalDist}{Logical. Whether the number of positions returned should be equally spaced or not}
}
\details{The distances between the new vessel positions are assumed to be equal.
}
\value{
A data frame in tacsat format is returned with estimated locations, speeds and headings between each real or actual vessel position.
The column added HL_ID refers to the ID given to each individual interpolation.
}
\references{EU Lot 2 project}
\author{Niels Hintzen, Doug Beare
}

\seealso{
\code{\link{interpolateTacsat}}
}
\examples{

data(tacsat)
# Take out first 100 rows
input   <- tacsat[1:100,]
# set resolution for number of locations or pings required
#   between each real or actual vessel location.
res1    <- 100

#Do the interpolation
interpolation <- interpolateTacsat(tacsat=input,interval=120,
                                   margin=12,res=res1,method="cHs",
                                   params=list(fm=0.4,distscale=20,sigline=0.2,
                                   st=c(4,8)),headingAdjustment=0)
                                          
#Convert the interpolated data back to tacsat with 10 locations or
# pings between each real or actual vessel location.
tacsat        <- interpolation2Tacsat(interpolation=interpolation,
                                      tacsat=input,npoints=10,equalDist=TRUE)
  
}
