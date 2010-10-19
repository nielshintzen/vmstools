\name{SplitLandingsAmongPings}
\alias{SplitLandingsAmongPings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distributes catches or values of catches among VMS pings

}
\description{
This function uses the output from countPings.r to distribute the catches of ONE chosen species recorded in eflalo over the pings in tacsat.
It then adds on the total weight and total landing values of all species in eflalo2.
\usage{
SplitLandingsAmongPings(tacsat = tacsat, eflalo2 = eflalo2, which.species = "COD")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tacsat}{
%%     data.frame, tacsat format 
}
  \item{eflalo2}{
%%     data.frame, eflalo format
}
  \item{which.species}{
%%     character, species name 
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
function(tacsat=tacsat, eflalo2=eflalo2, which.species="COD")
{
#This function uses the output from count pings to distribute the catches of ONE chosen species recorded in eflalo over the pings in tacsat.
#It then adds on the total weight and total landing values of all species in eflalo2
dn <- dimnames(eflalo2)[[2]]
yp <- NULL
for(ss in which.species){yp <- c(yp,grep(ss,dn)) }
kgs <- eflalo2[,yp[1]][match(tacsat$FT_REF,eflalo2$FT_REF)]
tacsat$kgs.by.ping <- kgs/tacsat$npings
cash <- eflalo2[,yp[2]][match(tacsat$FT_REF,eflalo2$FT_REF)]
tacsat$cash.by.ping <- cash/tacsat$npings
#All species
#KGS
yp <- grep("_KG_",dn)
eflalo2$kgs <- apply(eflalo2[,yp],1,sum,na.rm=T)
tacsat$total.kgs <- eflalo2$kgs[ match(tacsat$FT_REF,eflalo2$FT_REF)]
#VALUES
yp <- grep("_EURO_",dn)
eflalo2$cash <- apply(eflalo2[,yp],1,sum,na.rm=T)
tacsat$total.cash <- eflalo2$cash[ match(tacsat$FT_REF,eflalo2$FT_REF)]
tacsat

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
