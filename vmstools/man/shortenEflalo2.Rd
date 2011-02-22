\name{shortenEflalo2}
\alias{shortenEflalo2}
\title{Cuts up eflalo according to any particular combination of species to make it more manageable
}
\description{
Cuts up eflalo according to any particular combination of species to make it more manageable
}
\usage{
shortenEflalo2(data = eflalo, which.species = c("SPEC"))
}
\arguments{
  \item{data}{
eflalo formatted data
}
  \item{which.species}{
array of species names in FAO species names
}
}
\references{
EU lot 2 project
}
\author{
Doug Beare
}
\examples{
# The function is currently defined as
function(data=eflalo2,which.species=c("ANE","BIB","BLL","COD","DAB","HAD",
                                      "HER","MAC","NEP","PLE","SOL","WHG"))
{
#Eflalo is an unwieldy format so this is handy if you
#  want to select a few species
dn <- dimnames(data)[[2]]
yp <- NULL
for(ss in which.species){yp <- c(yp,grep(ss,dn)) }
short.eflalo2 <- data[,c(1:26,yp)]
short.eflalo2
  }
}