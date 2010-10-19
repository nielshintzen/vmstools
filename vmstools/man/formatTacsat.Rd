\name{formatTacsat}
\alias{formatTacsat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Format tacsat data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
formatTacsat(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function(x){
  x$VE_REF    <- ac(x$VE_REF)
  x$SI_LATI   <- an(x$SI_LATI)
  x$SI_LONG   <- an(x$SI_LONG)
  x$SI_DATE   <- ac(x$SI_DATE)
  x$SI_TIME   <- ac(x$SI_TIME)
  x$SI_SP     <- an(x$SI_SP)
  x$SI_HE     <- an(x$SI_HE)
  #Get rid of NAs in the long and lats
  x <- x[!is.na(x$SI_LATI),]
  x <- x[!is.na(x$SI_LONG),]
  return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
