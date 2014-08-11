\name{estimatePropFishing}
\alias{estimatePropFishing}
\title{
Estimate the proportion fishing per fishing trip
}
\description{
Estimate what proportion of time a fishing vessel is fishing when outside of harbour
}
\usage{
estimatePropFishing(tacsat, eflalo, by = c("LE_GEAR", "VE_REF"))
}
\arguments{
  \item{tacsat}{
tacsat dataset, including FT_REF column
}
  \item{eflalo}{
eflalo dataset
}
  \item{by}{
column names to calculate the estimation over
(e.g. LE_GEAR to estimate proportion fishing by gear or VE_REF to estimate
proportion fishing by vessel ID)
}
}
\value{
Returns the eflalo dataset including a new column 'propFish' which provides the
estimated proportion fishing for each entry}
\author{
Niels T. Hintzen}

\seealso{
\code{\link{raiseTacsat}}
}
