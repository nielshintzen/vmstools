\name{getMetiersClusters}
\alias{getMetiersClusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding Metiers from a reduced EFLALO dataset, step 3 : Clustering logevents using various multivariate methods 
}
\description{
This function represents the third step in the multivariate analysis of logbooks data for identifying metiers. 

This step allows applying various multivariate analyses on the data sets coming out of the first and second step. All methods will lead to 
a classification of all individuals (logevents), but they differ in their nature and then consequently in their outcomes. The four methods available are 
- Hierarchical Ascending Classification (HAC), with user-defined method for estimating distances and link for aggregating individuals 
- K-Means,
- Partitioning Around Medoids (PAM), 
- Clustering LARge Applications (CLARA).


The HAC method works by calculating the distance between individuals using the method selected with param1 ("euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski")
and aggregating them based on the distance between clusters, using the link selected with param2 ("ward", "single", "complete", "average", "mcquitty", "median" or "centroid").
In HAC, the number of classes is determined afterwards, once all combinations have been calculated, by using the objective criteria of scree test which
calculates the largest marginal loss of inertia between two consecutive numbers of classes. Therefore, the computing time and memory request for this method
can be quite comprehensive, and may reach memory limits on standard PC when operating with very large datasets. 

The K-Means method works by randomly choosing k individuals, or kernels (k corresponding to the final number of classes), and then affecting each individuals of the 
dataset to the closest kernel. Eachtime, the gravity center of the class is recalculated, thus reinitialising the calculation of distances to the next 
individual. In order to define the most appropriate number of classes, this procedure is repeated with differents values for k, from 2 to 15. The final number
of classes is identified by measuring the largest marginal loss of inertia between two consecutive numbers of classes. 

The PAM method works slightly around the same principle, starting with the initialisation of k medoids. The medoid is the individual in a class which shows
least dissimilarity with other individuals in the same class, and the remaining individuals are affected to their closest medoid. Then the sum of dissimilarities
is calculated and compared with the sum of dissimilarities if any other individual in the class had been playing the role of the medoid, and then the medoid
is eventually adjusted accordingly, until full stabilisation of the procedure. As in K-Means, the procedure is run with k varying from 2 to 15. The most appropriate
number of classes is identified by using the estimated silhouette of the classification for each value of k. The silhouette represent an average comparison of 
the distance between an individual and the other individuals from its class, and between the same individual and the other individuals from the next closest
class, and is therefore an objective measurement of the quality of the classification. The final number of classes retained is thus the one for which the silhouette
value is closest from 1. It is to be noted that the PAM method is not designed for working with very large datasets (>10 000 lines), and may quickly reach memory limits. 
 
The CLARA method is an extension of the PAM algorithm aiming at working with large datasets, and is therefore more efficient and powerful than PAM. 
It works by sampling 10 subsets of the dataset and running the PAM algorithm explained above on these subsets, and then to keep only the subset giving 
the best classification (silhouette). Afterwards, all remaining individuals are affected to their closest medoid, using user-defined method (param1) for calculating
distances ('euclidean' being used as the default, but 'manhattan' could also be used). 


}

\usage{
getMetierClusters(datSpecies,datLog,analysisName="",methMetier="clara",param1="euclidean",param2=NULL)

}

\arguments{
  \item{datSpecies}{numerical matrix with Logevents as lines and species as columns, with percentage values (between 0 and 100) of each species in the logevent catches. 
Logevent ID (LE_ID) should be as row names. Typically, this table will be produced from a eflalo dataset at the step 1 of the metier analysis, 
using the function extractTableMainSpecies()
}
  \item{datLog}{numerical matrix with Logevents as lines, and values to be used for calculating distances between individuals as columns. 
Typically, this table will be produced at the step 2 of the metier analysis, using the function getTableAfterPCA(). If a PCA was run, 
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
\details{
}

\value{
For any of the four methods used, a number of graphs are produced and saved directly in the working directory. 
They describe 1) projections of results on factorial plans, 2) the mean profile of the whole dataset, i.e. the average percentage of the various species in a 
logevent across all individuals, 3) the mean and standard deviation profile in terms of average percentage of catch by species within each cluster, 
4) the number of logevents by clusters, and 5) the profile of target species (TargetValue) by cluster. 

Secondly, the method will return an excel worksheet with one page per cluster, tabulating the main results for future use.

Finally, the function returns a list with a number of results and diagnotics on the performance of the method:
  \item{LE_ID_clust}{a data frame with two columns, linking the initial ID name of the Logevent (LE_ID) with the cluster number where the ID has been allocated.   
}
  \item{clusters}{diagnostics of the clustering process. It may vary between the four methods. 
}
  \item{nameTarget}{A tabulated list ranking the most characteristic (species with a test-value > 1.98 by cluster) species in the clusters. 
}
  \item{betweenVarClassifOnTot}{Percentage of variance explained by the classification.  
}
  \item{nbClust}{Final number of clusters retained.  
}
  \item{summaryClusters}{Array documenting, for each cluster, the minimum, mean, maximum, as well as the 25%, 50% and 75% quantiles values of the
  percentage of catch by species for the individual logevents in the cluster.     
}
  \item{resval}{Matrix of test values by species and cluster. The test value measures for each species the difference between the average percentage
  of catches in the cluster compared to the average percentage of catch in the total dataset, thus large positive values (>1.98) will point out the most characteristic species in the clusters.     
}
  \item{target}{A tabulated list ranking the most important (in terms of percentage value) and characteristic (in terms of test values) species in the clusters.     
}

 
}

\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sébastien Demanèche, Stepanie Mahévas, Clara Ulrich}
\note{A number of libraries are initially called for the whole metier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(RODBC) }

\seealso{\code{extractTableMainSpecies(),getTableAfterPCA()}}

\examples{

  \dontrun{
               
  }

}

\keyword{}
