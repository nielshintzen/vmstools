\name{compareToOrdination}
\alias{compareToOrdination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding Métiers from a reduced EFLALO dataset, step 4 : compare the métiers found with multivariate classification to the métiers
found with ordination method.
}
\description{
This function permits to compare the métiers found with multivariate classification to the métiers found with ordination method
for each logevent.

The ordination method consists to allocate the logevent to the first species in catch, and this species defines the metiers (level 7)
of the logevent. At level 5, the metiers is defined by the first group of species in catch.

The multivariate classification defines métiers (level 7), and each métier is characterized by one or several species,
contrary to the métiers from ordination method, which are described by only one species (level 7) or one group of species (level 5).
The métier (level 5) from multivariate classification is determined by the groups of species to which belong the
characteristic species of each métier. This method permits thus to obtain mixed métiers, defined by several groups of species
(level 5), what is impossible with the ordination method.

This function compareToOrdination() shows the differences between the métiers from ordination method and
from multivariate classification, and shows the interest of this last to obtain exhaustive and accurate métiers.

}

\usage{
compareToOrdination(dat, Step2, clusters)

}

\arguments{
  \item{dat}{a data.frame reduced from an eflalo format. It should contain only the LE_ID (Logevent ID) variable as well as
  all species names in columns, with raw catch data. It is necessary to sort out potential error-prone lines (such as lines
  with only 0) prior to the analysis, and to replace NA values by 0.
}
  \item{Step2}{numerical matrix with logevents as lines, and values to be used for calculating distances between individuals
   as columns. Typically, this table will be produced at the step 2 of the métier analysis, using the function getTableAfterPCA().
   If a PCA was run, selected Principal Components will appear as columns. If no PCA was run, the matrix will be the same
   as datSpecies (produced at the step 1 of the métier analysis, using the function extractTableMainSpecies()), with 
   percentage values by species.
}
  \item{clusters}{the clustering vector. An integer vector of length n, the number of logevents, giving for each
  logevents the number of the cluster to which it belongs. This vector will be produced at the step 3 of the métier analysis,
  using the function getMetierClusters().
}
}
\details{
}

\value{
The function returns a list with a number of tables on the comparison of the two methods :
  \item{nbLogFirstSp}{An integer vector, which gives the number of logevents allocate to each métier (level 7) defined by
  ordination method.
}
  \item{compNbLogClustVsFirstSp}{A table, which gives the distribution of the logevents of each métier (level 7) defined by 
  multivariate classification in the métiers (level 7) defined by ordination method.
}
  \item{compL57vsL5}{A table, which gives the distribution of the logevents of each métier (level 5) defined by multivariate
  classification in the métiers (level 5) defined by ordination method. 
}
  \item{compL57vsL7}{A table, which gives the distribution of the logevents of each métier (level 5) defined by multivariate
  classification in the métiers (level 7) defined by multivariate classification.
}
  \item{compL5vsL7}{A table, which gives the distribution of the logevents of each métier (level 5) defined by ordination
  method in the métiers (level 7) defined by multivariate classification.
}
}

\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sebastien Demanèche, Stepanie Mahévas, Clara Ulrich}
\note{A number of libraries are initially called for the whole metier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(RODBC) }

\seealso{\code{extractTableMainSpecies(), getTableAfterPCA(), getMetierClusters()}}

\examples{

  \dontrun{

  }

}

\keyword{}
