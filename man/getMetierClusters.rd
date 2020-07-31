\name{getMetierClusters}
\alias{getMetierClusters}
\title{
Finding metiers from a reduced EFLALO dataset, step 3: clustering logevents using various multivariate methods 
}
\description{
This function represents the third step in the workflow processing logbooks data for identifying metiers. 

This step allows applying various clustering analyses on the data sets coming out of the first and second step. All methods will lead to 
a classification of all individuals (logevents), but they differ in their nature and then consequently in their outcomes. The four methods available are 
- Hierarchical Ascending Classification (HAC), with user-defined method for estimating distances and link for aggregating individuals\cr
- K-Means,\cr
- Partitioning Around Medoids (PAM),\cr
- Clustering LARge Applications (CLARA).


The HAC method works by calculating the distance between individuals using the method selected with param1 ("euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski")
and aggregating them based on the distance between clusters, using the link selected with param2 ("ward", "single", "complete", "average", "mcquitty", "median" or "centroid").
In HAC, the number of classes is determined afterwards, once all combinations have been calculated, by using the objective criteria of scree test which
detects the third marginal loss of inertia between two consecutive numbers of classes. Therefore, the computing time and memory request for this method
can be quite comprehensive, and may reach memory limits on standard PC when operating with very large datasets. 

The K-Means method works by randomly choosing k individuals, or kernels (k corresponding to the final number of classes), and then affecting each individuals of the 
dataset to the closest kernel. Each time, the gravity center of the class is recalculated, thus reinitialising the calculation of distances to the next 
individual. In order to define the most appropriate number of classes, this procedure is repeated with differents values for k, from 2 to 15. The final number
of classes is identified by using the criteria of scree test which detects the third marginal loss of inertia between two consecutive numbers of classes.

The PAM method works slightly around the same principle, starting with the initialisation of k medoids. The medoid is the individual in a class which shows
least dissimilarity with other individuals in the same class, and the remaining individuals are affected to their closest medoid. Then the sum of dissimilarities
is calculated and compared with the sum of dissimilarities if any other individual in the class had been playing the role of the medoid, and then the medoid
is eventually adjusted accordingly, until full stabilisation of the procedure. The most appropriate number of classes is identified by using the estimated silhouette
of the classification for each value of k. The silhouette represent an average comparison of the distance between an individual and the other individuals from its
class, and between the same individual and the other individuals from the next closest class, and is therefore an objective measurement of the quality of the classification.
The final number of classes retained is the one for which the second maximum of the silhouettes is reached. It is to be noted that the PAM method is not designed
for working with very large datasets (>10 000 lines), and may quickly reach memory limits.
 
The CLARA method is an extension of the PAM algorithm aiming at working with large datasets, and is therefore more efficient and powerful than PAM. 
It works by sampling 5 subsets of the dataset for each value of k and running the PAM algorithm explained above on these subsets, and then to keep only the subset
giving the best classification (silhouette). As in PAM, the final number of classes retained is the one for which the second maximum of the silhouettes is reached.
Afterwards, all remaining individuals are affected to their closest medoid, using user-defined method (param1) for calculating distances ('euclidean' being used
as the default, but 'manhattan' could also be used).


}

\usage{
getMetierClusters(datSpecies,datLog,analysisName="",methMetier="clara",
                  param1="euclidean",param2=NULL)

}

\arguments{
  \item{datSpecies}{numerical matrix with logevents as rows and species as columns, with percentage values (between 0 and 100) of each species in the logevent catches. 
Logevent ID (LE_ID) should be as row names. Typically, this input table will be produced from the step 1 of the metier analysis applied on the eflalo initial data, 
using the function extractTableMainSpecies()
}
  \item{datLog}{numerical matrix with logevents as rows, and values to be used for calculating distances between individuals as columns. 
Typically, this input table is produced by the step 2 of the metier analysis, using the function getTableAfterPCA(). If a PCA was run, 
selected Principal Components will appear as columns. If no PCA was run, the matrix will be the same as datSpecies, with percentage values by species. 
}
  \item{analysisName}{character, the name of the run. Used for the file name of the plots. 
}
  \item{methMetier}{character. The name of the clustering method to be used. Must be chosen between "hac", "kmeans", "pam" and "clara".   
}
  \item{param1}{character. Parameter used for chosing the method calculating distances between individuals, to be used in HAC and CLARA algorithms. 
  For HAC, it can be chosen between "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". For CLARA, between "euclidean" and "manhattan". 
  For PAM and K-means, this must be set to NULL.  
}
  \item{param2}{character. Parameter used for chosing the method calculating distances between clusters, to be used in HAC algorithm. 
  For HAC, it can be chosen between "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".  
  For PAM, CLARA and K-means, this must be set to NULL.  
}

}
\value{
For any of the four methods used, a number of graphs are produced and saved directly in the working directory. 
They describe 1) projections of results on factorial plans, 2) the mean profile of the whole dataset, i.e. the average percentage of the various species in a 
logevent across all individuals, 3) the mean and standard deviation profile in terms of average percentage of catch by species within each cluster, 
4) the number of logevents by clusters, and 5) the profile of test-values by cluster. 


Finally, the function returns a list with a number of results and diagnostics on the performance of the method:
  \item{LE_ID_clust}{a data frame with two columns, linking the initial ID name of the Logevent (LE_ID) with the cluster metier where the ID has been allocated.
}
  \item{clusters}{diagnostics of the clustering process. It may vary between the four methods. 
}
  \item{betweenVarClassifOnTot}{percentage of variance explained by the classification.  
}
  \item{nbClust}{final number of clusters retained.  
}
  \item{summaryClusters}{array documenting, for each cluster, the minimum, mean, maximum, as well as the 25\%, 50\% and 75\% quantiles values of the
  percentage of catch by species for the individual logevents in the cluster.     
}
  \item{testValues}{matrix of test-values by species and cluster. The test-value measures for each species the difference between the average percentage
  of catches in the cluster compared to the average percentage of catch in the total dataset, thus large positive values (>1.98) will point out the most 
  characteristic species in the clusters.     
}
  \item{testValuesSpecies}{a tabulated list ranking the most characteristic species in the clusters (ordering species with a test-value > 1.98 by decreasing 
  order of test-value).     
}
  \item{percLogevents}{a matrix giving the percentage of logevents catching each species by cluster.
}
  \item{descClusters}{a data frame giving some descriptive statistics for each cluster, like cluster size, number of species needed to have at least 50\% of the
  cluster's total catch, number of species with a test-value > 1.98 in the cluster, number of species caught in at least 50\% of the logevents, etc...
}
  \item{tabClusters}{a 3d-array giving a table summurazing for each cluster the most important species (in terms of catch),
  the associated test-value, and the percentage of logevents of the cluster catching these species.     
}
  \item{targetSpecies}{a list giving the target species by cluster.
}
 
}

\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER, France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)}
\note{A number of libraries are initially called for the whole metier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)}

\seealso{\code{\link{getEflaloMetierLevel7}}, \code{\link{selectMainSpecies}},
\code{\link{extractTableMainSpecies}}, \code{\link{getMetierClusters}},
\code{\link{getTableAfterPCA}}}

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
  Step2 <- getTableAfterPCA(Step1, analysisName, pcaYesNo="pca",
              criterion="70percents")

  row.names(Step1) <- rowNamesSave
  row.names(Step2) <- rowNamesSave
  
  # Define a metier for each logevent running the CLARA algorithm 
  Step3 <- getMetierClusters(Step1, Step2, analysisName,
              methMetier="clara", param1="euclidean", param2=NULL)
              
  }

}
