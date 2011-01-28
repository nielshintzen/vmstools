\name{predictMetier}
\alias{predictMetier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding métiers for a new reduced EFLALO dataset from a reduced EFLALO dataset for which the métiers are known.
}
\description{
This function permits to determine the métiers of new logevents belonging to a new reduced EFLALO dataset, from a reduced EFLALO dataset for which the 
métiers are known. For example, this function permits to determine the métiers of logevents realized in 2008, basing to métiers of logevents realized 
in 2007. It works with a discriminante analysis ran with the function fda. In the first time, the function fda calibrates a model basing on the dataset 
for which the métiers are known. In the second time, the function predict is used to predict the métiers of new logevents, thanks to the previous 
model.
}

\usage{
predictMetier(learningData,clustersAffectation,newData)
}

\arguments{
  \item{learningData}{A numerical matrix with logevents as lines and species as columns, with percentage values (between 0 and 100) of each species 
  in the logevent catches. Logevent ID (LE_ID) should be as row names. Typically, this table will be produced from a eflalo dataset at the step 1 of 
  the métier analysis, using the function extractTableMainSpecies()
}
  \item{clustersAffectation}{An integer vector of length n, the number of logevents in learningData, giving for each
  logevent the number of the cluster to which it belongs. This vector will be produced at the step 3 of the métier analysis,
  using the function getMetierClusters(). 
}
  \item{newData}{A numerical matrix with logevents as lines and species as columns, with values of each species in the logevent catches. 
  Logevent ID (LE_ID) should be as row names. 
}
}

\details{
}

\value{
A data frame giving the logevent ID (LE_ID) (first column) and the number of the cluster to which it belongs (second column), for all logevents
belonging to newData.
}

\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sébastien Demanèche, Stéphanie Mahévas, Clara Ulrich}
\note{A number of libraries are initially called for the whole métier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)}

\seealso{\code{extractTableMainSpecies(),getMetierClusters()}}

\examples{

  \dontrun{
               
  }

}

\keyword{}
                                       