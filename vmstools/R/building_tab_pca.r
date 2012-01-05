#########################################
# Build the table with the main species #
#########################################

building_tab_pca=function(data,especes){
  p=ncol(data)
  noms=colnames(data)
  ind_princ=which(is.element(noms,especes))
  princ=data[,ind_princ]
  return(princ)
}
