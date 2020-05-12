\name{getTableAfterPCA}
\alias{getTableAfterPCA}
\title{
Finding Metiers from a reduced EFLALO dataset, step 2: Options for running a PCA on the selected species. 
}
\description{
This function represents the second step in the multivariate analysis of logbooks data for identifying metiers. 

This step of the classification is a user-driven choice of running or not a Principal Component Analysis (PCA)
on the catch per logevents of the species retained, expressed in percentage. A PCA will help reducing 
the multi-dimensional catch matrix to a smaller number of informative components represented by the first 
n- axes of the PCA transformation. The number of axes can be chosen either by using a scree test (criterion
'screetest') looking for the significant marginal increases of explained inertia, or by selecting all axes
cumulating 70 percent of explained inertia (criterion '70percents').

}

\usage{
getTableAfterPCA(datSpecies,analysisName="",
                 pcaYesNo="pca",criterion="70percents")

}

\arguments{
  \item{datSpecies}{numerical matrix with Logevents as lines and species as columns, with percentage values (between 0 and 100) of each species in the logevent catches. 
Logevent ID (LE_ID) should be as row names. Typically, this table will be produced from a eflalo dataset using the function extractTableMainSpecies()
}
  \item{analysisName}{character, the name of the run. Used for the file name of the plots. 
}
  \item{pcaYesNo}{character. An indication of whether a PCA should be performed. Use "pca" if a PCA should be run, and "nopca" if it shouldn't.  
}
  \item{criterion}{character. Criterion used for selected the number of principal components (axes) retained. 
  The number of axes can be chosen either by using a scree test (criterion 'screetest'), looking for the second-order unsignificant marginal increases 
  of explained inertia, or by selecting all axes cumulating up to 70 percent of explained inertia (criterion '70percents'). If pcaYesNo="nopca", then this criterion
  should be set as NULL. Experience has shown that the '70percents' criterion may often select a larger number axes than the scree test criterion, and may therefore be
  more appropriate for analysing large and heterogeneous datasets, by retaining more information for the subsequent clustering.  
}
}
\value{
If a PCA is run, a number of graphs are produced and directly saved in the working directory. These graphs describe 1) species projection on the first 
factorial axes, 2) individuals (logevents) projections on the first factorial axes, 3) actual and cumulative percentage of inertia of factorial axes, and 4) eigenvalues.  

If a PCA is run, the function returns a numerical matrix with Logevents as lines and selected Principal Components as columns. 
If no PCA is run, the function returns the same matrix as the input, with percentage values by species. 
}

\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER, France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)}
\note{A number of libraries are initially called for the whole metier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)}

\seealso{selectMainSpecies(), extractTableMainSpecies()}

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
  
  }

}
