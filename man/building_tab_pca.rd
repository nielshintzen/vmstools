\name{building_tab_pca}
\alias{building_tab_pca}


\title{
Useful functions for the multivariate analysis of logbooks data for identifying metiers.
}

\description{
This function contains several functions needed for the multivariate analysis of logbooks data for identifying metiers.
}
\arguments{
\item{transformation_proportion}{Transform quantities to percentage values (between 0 and 100) of each species in the logevent total catch.
}

\item{table_variables}{Transpose the dataset (change variables into individuals)
}

\item{scree}{Implementation of "scree-test"
}

\item{select_species}{Remove the cluster with the smallest mean of capture
}

\item{building_tab_pca}{Build the table with the main species
}

\item{test.values}{Compute the test-value for each species by cluster
}

\item{targetspecies}{Determine the species with a test-value > 1.96 by cluster
}

\item{withinVar}{Calculate the cluster's within-variance
}

}


\references{Development of tools for logbook and VMS data analysis. Studies for carrying out the common fisheries policy No MARE/2008/10 Lot 2}
\author{Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER, France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)}
\note{A number of libraries are initially called for the whole metier analyses and must be installed : (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)}
