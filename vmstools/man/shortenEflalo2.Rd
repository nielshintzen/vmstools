\name{shortenEflalo2}
\alias{shortenEflalo2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cuts up eflalo according to any particular combination of species to make it more manageable
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
shortenEflalo2(data = eflalo2, which.species = c("ANE", "BIB", "BLL", "COD", "DAB", "HAD", "HER", "MAC", "NEP", "PLE", "SOL", "WHG"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{which.species}{
%%     ~~Describe \code{which.species} here~~
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
Doug Beare
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
function(data=eflalo2,which.species = c("ANE","BIB","BLL","COD","DAB","HAD","HER","MAC","NEP","PLE","SOL","WHG") )
{
#Eflalo is an unwieldy format so this is handy if you want to select a few species
dn <- dimnames(data)[[2]]
yp <- NULL
for(ss in which.species){yp <- c(yp,grep(ss,dn)) }
short.eflalo2 <- data[,c(1:26,yp)]
short.eflalo2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
