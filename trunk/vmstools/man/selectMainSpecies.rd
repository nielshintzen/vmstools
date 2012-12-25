\name{selectMainSpecies}
\alias{selectMainSpecies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identyfing in an EFLALO dataset which species can be considered as important for the analysis of target species, by crossing three different approaches. 
}
\description{
A typical logbook dataset will contain a large number of species recorded, but only a limited 
number of these could potentially be considered as target species. This function aims thus at identifying these
by using three different approaches :
   - HAC (Hierarchical Ascending Classification) based on Euclidian distances between species with Ward aggregating criteria;
   - Total, where species are ranked based on their proportion in the total catch, and those cumulating to a given percentage are retained
   - Logevent, where species are selected if they represent at least a given percentage of the logevent's catch for at least one logevent (one line)

The HAC and the Logevent methods work on catch data transformed in percentage of species by logevent (line), in order to remove the effect of large
hauls compared to small hauls. The Total method works on raw data.     
In the HAC method, a first group of species, the residual ones, is identified by clustering and using a first-order scree test for cutting the tree. 
Other species are pooled in the group of principals. New similar HACs are run through a loop on this group of residuals species, to identify if any new species
might have been left aside in the first run. It is important to note though that HAC method might quickly reach memory limits on standard PCs, and may thus
not be run on very large datasets.
In the Total method, the percentage threshold is being increased with 5\% steps from 5 to 100, and the ranked species summing up to this value
is recorded. 
In the Logevent method, the percentage threshold is also being increased with 5\% steps from 5 to 100, and all species representing at least this value in
at least one line are recorded. 

This function allows thus to explore the variability and the sensitivity of the definition of key species to differences in concepts and subjective
thresholds. A plot showing the number of species retained according to these two criteria is produced, allowing the user to make a qualitative choice
when this will be used for identyfing metiers. Empirical experience led us to suggest that the combination of species entering either the HAC, or 95\% of the
total cumulated catch, or 100\% of at least one logevent would be sufficient to cover the main target species of a standard logbook dataset, but other
choices could be made.          

  

}

\usage{
selectMainSpecies(dat,analysisName="",RunHAC=TRUE,DiagFlag=FALSE)

}

\arguments{
  \item{dat}{a data.frame reduced from an eflalo format. It should contain only the LE_ID (Logevent ID) variable as well as all species names in columns, 
  with raw catch data. It is necessary to sort out potential error-prone lines (such as lines with only 0) prior to the analysis, and to replace NA values by 0.
}
  \item{analysisName}{character, the name of the run. Used for the file name of the plots. 
}
  \item{RunHAC}{Boolean. In case of very large datasets, memory limits might be reached with HAC method. This option allows thus to skip the method 
  and consider species selection using the two other methods only. 
}
  \item{DiagFlag}{Boolean. If DiagFlag=TRUE, additional plots and diagnostics are produced. Not very used. 
}
}
\value{
The function produces a plot (saved in the "analysisName" working directory showing the number of species selected according to the method and 
the percentage threshold selected for both 'Total' and 'Logevent' methods. 
The function returns also a list of diagnostics of the three methods :
  \item{nbAllSpecies}{Number of species initially in the dataset 
}
  \item{propNbMainSpeciesHAC}{Proportion of the number of species retained by the HAC method to the total number of species 
}
  \item{propNbMainSpeciesTotal}{Proportion of the number of species retained by the Total method to the total number of species 
}
  \item{propNbMainSpeciesLogevent}{Proportion of the number of species retained by the Logevent method to the total number of species 
}
  \item{nbMainSpeciesHAC}{Number of species retained by the HAC method  
}
  \item{nbMainSpeciesTotal}{Number of species retained by the Total method with a percentage threshold increasing at 5\% step from 5\% to 100\%
}
  \item{nbMainSpeciesLogevent}{Number of species retained by the Logevent method with a percentage threshold increasing at 5\% step from 5\% to 100\%
}
  \item{namesMainSpeciesHAC}{Names of species retained by the HAC method 
}
  \item{namesMainSpeciesTotalAlphabetical}{Names of species retained by the Total method with 95\% threshold in alphabetical order, for easier comparison
  with the two other methods
}
  \item{namesMainSpeciesTotalByImportance}{Names of species retained by the Total method with 95\% threshold in ranked order of importance
}
  \item{namesMainSpeciesLogevent}{Names of species retained by the Logevent method with 100\% threshold
}
  \item{namesMainSpeciesAll}{Unique combination of the species retained in either HAC method, Total method with 95\% threshold, and Logevent method
  with 100\% threshold
}
  \item{propCatchMainSpeciesAll}{Proportion of the total catch represented by the selected species (species in namesMainSpeciesAll)
}
}

\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER, France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)}
\note{A number of libraries are initially called for the whole metier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)}

\seealso{\code{\link{extractTableMainSpecies}}}

\examples{


data(eflalo)

eflalo <- formatEflalo(eflalo)

eflalo <- eflalo[eflalo$LE_GEAR=="OTB",]

# note that output plots will be sent to getwd()
analysisName <- "metier_analysis_OTB"

dat <- eflalo[,c("LE_ID",grep("EURO",colnames(eflalo),value=TRUE))]
names(dat)[-1] <- unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))

explo <- selectMainSpecies(dat, analysisName, RunHAC=TRUE, DiagFlag=FALSE)
  #=> send the LE_ID and LE_EURO_SP columns only
}
