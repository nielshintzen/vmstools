\name{countPings}
\alias{countPings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ping counting over trips
}
\description{
This function counts the number of pings in the tacsat data for each trip and then merges them with that same file.
}
\usage{
countPings(tacsat = tacsat1f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tacsat}{
%%     data.frame, tacsat typically with only fishing pings identified
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(tacsat=tacsat1f)
{
#This function counts the number of pings in the tacsat data for each trip and then merges them with that same file.
npings.by.trip <-  aggregate(tacsat$FT_REF,list(tacsat$FT_REF),length)
dimnames(npings.by.trip)[[2]] <- c("FT_REF","npings")
tacsat.with.pings <- merge(tacsat,npings.by.trip,all=T)
tacsat.with.pings
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
