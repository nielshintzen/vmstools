########################################################
# Remove the cluster with the smallest mean of capture #
########################################################

select_species=function(data,groupes_cah){
  nb.classes=length(levels(as.factor(groupes_cah)))
  moyennes=numeric(nb.classes)
  for(i in 1:nb.classes){
    namegp=names(which(groupes_cah==i))
    effgp=length(which(groupes_cah==i))
    moyennes[i]=sum(data[namegp],na.rm=TRUE)/effgp
  }
  indice.autre=which(moyennes == min(moyennes,na.rm=TRUE))
  noms=names(which(groupes_cah!=indice.autre))
  return(list(noms,indice.autre))
}

