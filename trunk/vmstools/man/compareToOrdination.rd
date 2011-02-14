\name{compareToOrdination}
\alias{compareToOrdination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare the métiers found with the classification from clustering e.g. CLARA
against the métiers found by two simple alternative ordination methods: 'first species' &  'first group'.
}
\description{
This function allows to compare the métiers found using the classification from clustering 
against the métiers found by two alternative ordination methods, 'first species' and 'first group'.

The 'first species' method consists of characterizing the level 7 of each logevent by the species with the highest catch. 
Then, the level 5 is deduced from the target assemblage of species for which this species belongs to.

The 'first group' method alternatively assigns a level 5 to each logevent correponding to the target assemblage having the 
highest catch volume.

The classification from clustering defines level 7 métiers for which each metier is characterized by one or several species,
 called target species, unlike the métiers from the simpler ordination methods. Hence, mixed métiers can be obtained 
 in case that these species belong to different target assemblages.

This function enlights the differences when defining the métiers between the classification from clustering and the simpler
ordination methods, and potentially demonstrates the higher power of the data clustering method in obtaining 
exhaustive and accurately defined métiers.

}

\usage{
compareToOrdination(dat, Step2, clusters, tabClusters)

}

\arguments{
  \item{dat}{a reduced data.frame from an eflalo format. It should contain only the LE_ID (Logevent ID) variable as well as
  all species names in columns, with raw catch data. It may be needed to sort out potential error-prone rows (such as rows
  with only 0) prior to the analysis, and to also replace NA values by 0.
}
  \item{Step2}{numerical matrix with logevents as rows, and values to be used for calculating distances between individuals
  as columns. This matrix is produced at the step 2 of the métier analysis, output of the function getTableAfterPCA().
  In case a PCA is run, the selected axes will appear as columns. If no PCA is run, the matrix will be the same
  as datSpecies (produced at the step 1 of the métier analysis, using the function extractTableMainSpecies()), with 
  percentage values by species.
}
  \item{clusters}{the vector storing the cluster label of each logbook event. This vector will be produced at the step 3 of the métier analysis,
  using the function getMetierClusters().
}
  \item{tabClusters}{a 3d-array giving a table summurazing for each cluster the most important species (in terms of catch),
  the associated test-value, and the percentage of logevents of the cluster catching these species.  
}
}
\details{
}

\value{
The function returns a list with a number of tables on the comparison of the three methods:
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
\author{Nicolas Deporte, Sébastien Demanèche, Stéphanie Mahévas, Clara Ulrich, Francois Bastardie}
\note{A number of libraries are initially called for the whole métier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda) }

\seealso{\code{getEflaloMetierLevel7(), selectMainSpecies(), extractTableMainSpecies(), getMetierClusters(), getTableAfterPCA(), getMetierClusters()}}

\examples{

  
  
  \dontrun{
  
  data(eflalo)
  
  eflalo <- formatEflalo(eflalo)

  eflalo <- eflalo[eflalo$LE_GEAR=="OTB",]

  analysisName <- "metier_analysis_OTB" # note that output plots will be sent to getwd()
  
  explo <- selectMainSpecies(
             dat=eflalo[,c("LE_ID",grep("EURO",colnames(eflalo),value=T))],
               analysisName, RunHAC=TRUE, DiagFlag=FALSE)
    #=> send the LE_ID and LE_KG_SP columns only
           
  Step1 <- extractTableMainSpecies(
              eflalo[,c("LE_ID",grep("EURO",colnames(eflalo),value=T))],
                 explo$NamesMainSpeciesHAC, paramTotal=95, paramLogevent=100)
    #=> send the LE_ID and LE_KG_SP columns only             

  # Run a PCA
  Step2 <- getTableAfterPCA(Step1, analysisName, pcaYesNo="pca", criterion="70percents")

  # Define a metier for each logevent running the CLARA algorithm 
  Step3 <- getMetierClusters(Step1, Step2, analysisName, methMetier="clara", param1="euclidean", param2=NULL)
  
  # Compare the differences between the métiers defined by CLARA (data clustering)
  # and the métiers defined by two simple ordination methods
  compMetiers <- compareToOrdination(
                   dat=eflalo[,c("LE_ID",grep("EURO",colnames(eflalo),value=T))],
                           Step2=Step2, clusters=Step3$clusters$clustering, tabClusters=Step3$tabClusters)
    #=> send the LE_ID and LE_KG_SP columns only             

  }
  
}

\keyword{}
