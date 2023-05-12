######################################################
# Compute the test-value for each species by cluster #
######################################################



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
#' @export test.values
test.values=function(groupes,data){

  n=nrow(data)
  p=ncol(data)
  noms_var=colnames(data)
  nb_groupes=length(levels(as.factor(groupes)))
  noms_groupes=character(nb_groupes)

  stats_globales=matrix(0,nrow=p,ncol=2)
  row.names(stats_globales)=noms_var
  colnames(stats_globales)=c("mean","variance")
  for (i in 1:p){
    stats_globales[i,1]=mean(data[,noms_var[i]])
    stats_globales[i,2]=var(data[,noms_var[i]])
  }

  res=matrix(0,nrow=p,ncol=nb_groupes)
  row.names(res)=noms_var

  for (j in 1:nb_groupes){
    groupe=which(groupes==j)
    n_k=length(groupe)

    for (i in 1:p){
      mu_k=mean(data[groupe,noms_var[i]])
      mu=stats_globales[noms_var[i],"mean"]
      V=stats_globales[noms_var[i],"variance"]
      V_mu_k=(n-n_k)*V/(n_k*(n-1))

      if(V_mu_k==0){
        Valeur_test=0
      }else{
        Valeur_test=(mu_k-mu)/sqrt(V_mu_k)
      }

      res[i,j]=Valeur_test
      rm(Valeur_test)
    }
    rm(groupe)
    noms_groupes[j]=paste("Cluster",j,sep=" ")
  }
  colnames(res)=noms_groupes
  return(res)
}
