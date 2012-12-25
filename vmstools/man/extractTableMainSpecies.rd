\name{extractTableMainSpecies}
\alias{extractTableMainSpecies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding Metiers from a reduced EFLALO dataset, step 1 : selecting the species considered as important for the analysis of target assemblages, 
by crossing three different approaches. 
}
\description{
This function represents the first step in the multivariate analysis of logbooks data for identifying metiers. 

A typical logbook dataset will contain a large number of species recorded, but only a limited 
number of these could potentially be considered as target species. This function builds further on the exploration conducted with the 
function selectMainSpecies() by combining three different approaches :
   - HAC (Hierarchical Ascending Classification) based on Euclidian distances between species with Ward aggregating criteria;
   - Total, where species are ranked based on their proportion in the total catches, and those cumulating to a given percentage are retained
   - Logevent, where species are selected if they represent at least a given percentage of at least one logevent (one line)

All the three methods work on catch data transformed in percentage of species by logevent (line), in order to remove the effect of large
hauls compared to small hauls.    
The HAC method is not run in this function, but in the selectMainSpecies() function, which outputs are directly used here. In addition, 
the relevant thresholds for the two other methods Total and Logevent must be selected.  

}

\usage{
extractTableMainSpecies(dat,namesMainSpeciesHAC,
                        paramTotal=95,paramLogevent=100)

}

\arguments{
  \item{dat}{a data frame reduced from an eflalo format. It should contain only the LE_ID (Logevent ID) variable as well as all species names in columns, 
  with raw catch data. It is necessary to sort out potential error-prone lines (such as lines with only 0) prior to the analysis, and to replace NA values by 0.
}
  \item{namesMainSpeciesHAC}{character. The names of species retained by the HAC method selectMainSpecies() function. If omitted, only the two other methods will be used for selecting species. 
}
  \item{paramTotal}{numeric. The percentage threshold for the Total method. All species ranked in terms of percentage in total catches and cumulating up to this 
  percentage value will be retained. Default value is 95\%.
}
  \item{paramLogevent}{numeric. The percentage threshold for the Logevent method. All species representing at least this threshold of at least one logevent (one line)
  will be retained. Default value is 100\% (logevents catching one and only one species).
}
}
\value{
The function returns a numerical matrix with Logevents as lines and species as columns, with percentage values (between 0 and 100) of each species in the logevent catches. 
Logevent ID (LE_ID) are returned as row names. 
}

\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER, France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)}
\note{A number of libraries are initially called for the whole metier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)}

\seealso{\code{\link{selectMainSpecies}}}

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
         
Step1 <- extractTableMainSpecies(dat, explo$namesMainSpeciesHAC,
            paramTotal=95, paramLogevent=100)
#=> send the LE_ID and LE_EURO_SP columns only
             
}

