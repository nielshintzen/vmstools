\name{pings2Fishframe}
\alias{pings2Fishframe}
\title{convert the merged (VMS + logbooks) data.frame to comply with the fishframe format}
\description{
the merged data.frame coming from mergeEflalo2Pings()
can be further converted to be uploaded in the fishframe data warehouse.
This includes the aggregation of fishing pings (VMS effort and logbook landings) 
by CSquare area code and month. 
}
\usage{
pings2Fishframe (general=list(output.path=file.path("C:","output"), a.year=2009, a.country="NLD") )}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{general}{a list with general settings, e.g. the output.path}
}
\details{
The function first starts by searching the input data.frame in the output.path folder.
The result will be saved in the output.path folder as well.
the 'data.table' package is required to perform the aggregation of the
data.frame much faster.
An area code is also added (e.g. from the ICESarea code)
}
\source{
}
\references{
}
\author{
Francois Bastardie
}

\examples{
  
  data(tacsat)
  
  # CONVERT TO FISHFRAME FORMAT (might take some time running)
  # (by default, this will keep all the species in the output table)
  tmp <- bindAllMergedTables (vessels= unique(tacsat$VE_REF), species.to.keep=character(), 
                      folder = file.path("C:","output"), all.in.one.table=FALSE)
  
  ff  <- pings2Fishframe (general=list(output.path=file.path("C:","output"),
                                                   a.year=2009, a.country="NLD") )

                                        

}
