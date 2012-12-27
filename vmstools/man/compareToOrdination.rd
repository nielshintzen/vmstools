\name{compareToOrdination}
\alias{compareToOrdination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare the metiers found with the classification from clustering e.g. CLARA
against the metiers found by two simple alternative ordination methods: 'first species' &  'first group'.
}
\description{
This function allows to compare the metiers found using the classification from clustering
against the metiers found by two alternative ordination methods, 'first species' and 'first group'.

The 'first species' method consists of characterizing the level 7 of each logevent by the species with the highest catch. 
Then, the level 5 is deduced from the target assemblage of species for which this species belongs to.

The 'first group' method alternatively assigns a level 5 to each logevent correponding to the target assemblage having the 
highest catch volume.

The classification from clustering defines level 7 metiers for which each metier is characterized by one or several species,
 called target species, unlike the metiers from the simpler ordination methods. Hence, mixed metiers can be obtained
 in case that these species belong to different target assemblages.

This function enlights the differences when defining the metiers between the classification from clustering and the simpler
ordination methods, and potentially demonstrates the higher power of the data clustering method in obtaining 
exhaustive and accurately defined metiers.

}

\usage{
compareToOrdination(dat, Step2, clusters, targetSpecies)

}

\arguments{
  \item{dat}{a reduced data.frame from an eflalo format. It should contain only the LE_ID (Logevent ID) variable as well as
  all species names in columns, with raw catch data. It may be needed to sort out potential error-prone rows (such as rows
  with only 0) prior to the analysis, and to also replace NA values by 0.
}
  \item{Step2}{numerical matrix with logevents as rows, and values to be used for calculating distances between individuals
  as columns. This matrix is produced at the step 2 of the metier analysis, output of the function getTableAfterPCA().
  In case a PCA is run, the selected axes will appear as columns. If no PCA is run, the matrix will be the same
  as datSpecies (produced at the step 1 of the metier analysis, using the function extractTableMainSpecies()), with
  percentage values by species.
}
  \item{clusters}{the vector storing the cluster label of each logbook event. This vector will be produced at the step 3 of the metier analysis,
  using the function getMetierClusters().
}
  \item{targetSpecies}{a list giving the target species by cluster. This list will be produced at the step 3 of the metier analysis,
  using the function getMetierClusters(). 
}
}

\value{
The function returns a list with a number of tables on the comparison of the three methods:
  \item{nbLogFirstSp}{An integer vector giving the number of logevents allocated to each metier (level 7) defined by
  the ordination method 'first species'.
}
  \item{compClustersL7vsFirstSpeciesL7}{A table giving the distribution of the logevents of each metier (level 7) defined by
  the multivariate classification in the metiers (level 7) defined by the ordination method 'first species'.
}
  \item{compClustersL5vsFirstSpeciesL5}{A table giving the distribution of the logevents of each metier (level 5) defined by
  the multivariate classification in the metiers (level 5) defined by the ordination method 'first species'.
}
  \item{compClustersL5vsFirstGroupL5}{A table giving the distribution of the logevents of each metier (level 5) defined by
  the multivariate classification in the metiers (level 5) defined by the ordination method 'first group'.
}
  \item{compFirstSpeciesL5vsFirstGroupL5}{A table giving the distribution of the logevents of each metier (level 5) defined by
  the ordination method 'first species' in the metiers (level 5) defined by the ordination method 'first group'.
}
  \item{compClustersL5vsClustersL7}{A table giving the distribution of the logevents of each metier (level 5) defined by
  the multivariate classification in the metiers (level 7) defined by the multivariate classification.
}
}

\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER, France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)}
\note{A number of libraries are initially called for the whole metier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)}

\seealso{\code{\link{getEflaloMetierLevel7}}, \code{\link{selectMainSpecies}},
\code{\link{extractTableMainSpecies}}, \code{\link{getMetierClusters}},
\code{\link{getTableAfterPCA}}, \code{\link{getMetierClusters}}}

\examples{

  
  
  \dontrun{
  
  data(eflalo)
  
  eflalo <- formatEflalo(eflalo)

  eflalo <- eflalo[eflalo$LE_GEAR=="OTB",]

  # note that output plots will be sent to getwd()
  analysisName <- "metier_analysis_OTB"
  
  dat <- eflalo[,c("LE_ID",grep("EURO",colnames(eflalo),value=TRUE))]
  names(dat)[-1] <- unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))

  explo <- selectMainSpecies(dat, analysisName, RunHAC=TRUE, DiagFlag=FALSE)
    #=> send the LE_ID and LE_KG_SP columns only
           
  Step1 <- extractTableMainSpecies(dat, explo$NamesMainSpeciesHAC,
              paramTotal=95, paramLogevent=100)
    #=> send the LE_ID and LE_KG_SP columns only             

  rowNamesSave <- row.names(Step1)
  row.names(Step1) <- 1:nrow(Step1)
  
  # Run a PCA
  Step2 <- getTableAfterPCA(Step1, analysisName,
              pcaYesNo="pca", criterion="70percents")

  row.names(Step1) <- rowNamesSave
  row.names(Step2) <- rowNamesSave

  # Define a metier for each logevent running the CLARA algorithm 
  Step3 <- getMetierClusters(Step1, Step2, analysisName,
              methMetier="clara", param1="euclidean", param2=NULL)
  
  # Compare the differences between the metiers defined by CLARA
  # and the metiers defined by two simple ordination methods
  compMetiers <- compareToOrdination(dat, Step2,
                    clusters=Step3$clusters$clustering,
                    targetSpecies=Step3$targetSpecies)
    #=> send the LE_ID and LE_KG_SP columns only             

  }
  
}
