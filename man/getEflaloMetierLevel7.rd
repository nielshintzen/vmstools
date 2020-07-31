\name{getEflaloMetierLevel7}
\alias{getEflaloMetierLevel7}
\title{
Finding Metiers from a full EFLALO dataset, general function calling all the various steps of the analyses at once
}
\description{
This function will run the various steps of the multivariate analyses in the right order, allowing the full anlysis to be run at once.


}

\usage{
getEflaloMetierLevel7(dat,analysisName,path,critData="EURO",
    runHACinSpeciesSelection=TRUE,paramTotal=95,paramLogevent=100,
    critPca="PCA_70",algoClust="CLARA")

}

\arguments{
  \item{dat}{a data.frame in eflalo format. It should be previously subsetted for the desired region and gear of interest.
}
  \item{analysisName}{character, the name of the run. Used for the file name of the plots.
}
  \item{path}{character, the path of the main directory where all data and plots will be stored. A number of sub-directories will be created by the function.
}
  \item{critData}{character. Choice if the analysis must be performed on eflalo data in weight ("KG") or in value ("EURO")
}
  \item{runHACinSpeciesSelection}{Boolean. the HAC can  be time consuming for large datasets, and can therefore be switched off for quicker runs.
}
  \item{paramTotal}{numeric. The percentage threshold for the Total method. All species ranked in terms of percentage in total catches and cumulating up to this
  percentage value will be retained. Default value is 95\%.
}
  \item{paramLogevent}{numeric. The percentage threshold for the Logevent method. All species representing at least this threshold of at least one logevent (one line)
  will be retained. Default value is 100\% (logevents catching one and only one species).
}
  \item{critPca}{character. Choice if a PCA must be run or not (="NO_PCA"), and if yes, which criterion is used for selecting the number of axes, to be chosen between
  scree test (="PCA_SC") or 70\% of cumulated inertia (="PCA_70")
}
  \item{algoClust}{character. Choice of the clustering algorithm,
                   either "CLARA", "HAC" or "KMEANS".
}
}

\value{

The function returns the same eflalo data set as was input, with an additional column giving the cluster name where each logevent has been allocated.

Beside, a number of sub-directories are created, with intermediate results and plots.



}

\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER, France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)}
\note{A number of libraries are initially called for the whole metier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)}

\seealso{
\code{\link{selectMainSpecies}},
\code{\link{extractTableMainSpecies}},
\code{\link{getTableAfterPCA}},
\code{\link{getMetierClusters}},
}
\examples{

\dontrun{

data(eflalo)

eflalo <- formatEflalo(eflalo)

eflalo <- eflalo[eflalo$LE_GEAR=="OTB",]

# Note that output plots will be sent to getwd()
path <- getwd()
analysisName <- "metier_analysis_OTB"

# Return a fully working eflalo dataset with a metier at DCF Level7 for each
# logbook event
eflalo_metiers_level7 <- getEflaloMetierLevel7(eflalo, analysisName, path,
    critData="EURO",runHACinSpeciesSelection=TRUE, paramTotal=95,
    paramLogevent=100,critPca="PCA_70", algoClust="CLARA")
        
}

}
