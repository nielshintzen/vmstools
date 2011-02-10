\name{compareToOrdination}
\alias{compareToOrdination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding Métiers from a reduced EFLALO dataset, step 4 : compare the métiers found with the multivariate classification to the métiers
found with two ordination methods : ordination method 'first species' & ordination method 'first group'.
}
\description{
This function permits to compare the métiers found with the multivariate classification to the métiers found with two ordination methods
for each logevent.

The ordination method 'first species' consists to allocate the logevent to the first species in catch, and this species defines the métiers (level 7)
of the logevent. Then, the group of species (level 5), to which the species belongs, defines the métier (level 5) of the logevent.

The ordination method 'first group' works directly at the level 5, and the métier (level 5) is defined by the first group of species in catch.

The multivariate classification defines métiers (level 7), and each métier is characterized by one or several species, called target species,
contrary to the métiers from ordination methods, which are described by only one species (level 7) ('first species') or one group of species (level 5)
('first group').
Then, the métier (level 5) from the multivariate classification is determined by the groups of species to which belong the
target species of each métier. This method permits thus to obtain mixed métiers, defined by several groups of species
(level 5), what is impossible with the ordination methods.

This function compareToOrdination() shows the differences between the métiers from the ordination methods and
from the multivariate classification, and shows the interest of the latter to obtain exhaustive and accurate métiers.

}

\usage{
compareToOrdination(dat, Step2, clusters, tabClusters)

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
  \item{tabClusters}{a three-dimensional array giving for each cluster a table resuming the most important (in catch) and 
  characteristic (in test-value) species in the cluster, and the percentage of logevents of the cluster catching these species.  
}
}
\details{
}

\value{
The function returns a list with a number of tables on the comparison of the three methods :
  \item{nbLogFirstSp}{An integer vector giving the number of logevents allocated to each métier (level 7) defined by
  the ordination method 'first species'.
}
  \item{compClustersL7vsFirstSpeciesL7}{A table giving the distribution of the logevents of each métier (level 7) defined by 
  the multivariate classification in the métiers (level 7) defined by the ordination method 'first species'.
}
  \item{compClustersL5vsFirstSpeciesL5}{A table giving the distribution of the logevents of each métier (level 5) defined by 
  the multivariate classification in the métiers (level 5) defined by the ordination method 'first species'. 
}
  \item{compClustersL5vsFirstGroupL5}{A table giving the distribution of the logevents of each métier (level 5) defined by 
  the multivariate classification in the métiers (level 5) defined by the ordination method 'first group'.
}
  \item{compFirstSpeciesL5vsFirstGroupL5}{A table giving the distribution of the logevents of each métier (level 5) defined by 
  the ordination method 'first species' in the métiers (level 5) defined by the ordination method 'first group'.
}
  \item{compClustersL5vsClustersL7}{A table giving the distribution of the logevents of each métier (level 5) defined by 
  the multivariate classification in the métiers (level 7) defined by the multivariate classification.
}
}

\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sébastien Demanèche, Stéphanie Mahévas, Clara Ulrich}
\note{A number of libraries are initially called for the whole métier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda) }

\seealso{\code{selectMainSpecies(), extractTableMainSpecies(), getTableAfterPCA(), getMetierClusters()}}

\examples{

  library(vmstools)
  
  # Load a dataset reduced from an eflalo format
  data(eflalo)
  analysisName="metiers_analysis"
  
  # Select the main species for the analysis
  explo=selectMainSpecies(eflalo,analysisName,RunHAC=TRUE,DiagFlag=FALSE)
  Step1=extractTableMainSpecies(eflalo,explo$NamesMainSpeciesHAC,paramTotal=95,paramLogevent=100)
  
  # Run a PCA
  Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="70percents")
  
  # Find the metier for each logevent thanks to the algorithm "CLARA"
  Step3=getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="clara",param1="euclidean",param2=NULL)
  
  # Compare the differences between the métiers defined by CLARA (multivariate classification) and the métiers defined by the two ordination methods
  compMetiers=compareToOrdination(dat=eflalo,Step2=Step2,clusters=Step3$clusters$clustering,tabClusters=Step3$tabClusters)
  
  
  \dontrun{

  }
  
}

\keyword{}
