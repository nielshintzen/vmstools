#############################################################
# Determine the species with a test-value > 1.96 by cluster #
#############################################################



#' Useful functions for the multivariate analysis of logbooks data for
#' identifying metiers.
#' 
#' This function contains several functions needed for the multivariate
#' analysis of logbooks data for identifying metiers.
#' 
#' 
#' @param transformation_proportion Transform quantities to percentage values
#' (between 0 and 100) of each species in the logevent total catch.
#' @param table_variables Transpose the dataset (change variables into
#' individuals)
#' @param scree Implementation of "scree-test"
#' @param select_species Remove the cluster with the smallest mean of capture
#' @param building_tab_pca Build the table with the main species
#' @param test.values Compute the test-value for each species by cluster
#' @param targetspecies Determine the species with a test-value > 1.96 by
#' cluster
#' @param withinVar Calculate the cluster's within-variance
#' @note A number of libraries are initially called for the whole metier
#' analyses and must be installed :
#' (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)
#' @author Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER,
#' France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)
#' @references Development of tools for logbook and VMS data analysis. Studies
#' for carrying out the common fisheries policy No MARE/2008/10 Lot 2
#' @export targetspecies
targetspecies=function(resval){
  p=nrow(resval)
  nbgp=ncol(resval)

  tabnumespcib=data.frame()
  tabnomespcib=data.frame()

  for(i in 1:nbgp){
    # qnorm(0.975,mean=0,sd=1)=1.96     (P(resval>1.96)=0.025)
    numespcib=which(resval[,i]>1.96)
    numespcibdec=numespcib[order(resval[numespcib,i],decreasing=TRUE)]
    nomespcib=names(numespcibdec)

    nbespgpcib=length(numespcib)

    if(nbespgpcib>0){
      for (j in 1:nbespgpcib){
        tabnumespcib[i,j]=numespcibdec[j]
        tabnomespcib[i,j]=nomespcib[j]
      }
    }else{
        tabnumespcib[i,]=NA
        tabnomespcib[i,]=NA
    }
  }
  tabnumespcib=as.matrix(tabnumespcib)
  tabnomespcib=as.matrix(tabnomespcib)
  return(list(tabnumespcib=tabnumespcib,tabnomespcib=tabnomespcib))
}

