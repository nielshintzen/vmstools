statCluster=function(cahacpclara) {

# Descriptives statistics about clusters
for(i in 1:cahacpclara$nbClust){
  print(paste("Cluster ",i))
  print(length(which(cumsum(cahacpclara$mprofil[i,][order(cahacpclara$mprofil[i,],decreasing=T)])<50))+1)
  print(length(which(cumsum(cahacpclara$mprofil[i,][order(cahacpclara$mprofil[i,],decreasing=T)])<90))+1)
  print(length(which(cahacpclara$resval[,i]>1.98)))
  print(length(which(cahacpclara$resval[,i]>3.29)))
#  print(length(which(apply(datSpecies,2,function (x) (sizeClusters[i]-length(which(x[clusters==i]==0)))/sizeClusters[i]*100)>50)))
#  print(length(which(apply(datSpecies,2,function (x) (sizeClusters[i]-length(which(x[clusters==i]==0)))/sizeClusters[i]*100)>90)))
#  print(length(which(apply(datSpecies,2,function (x) (effgroupe[i]-length(which(x[clusters==i]==0)))/effgroupe[i]*100)>50)))
#  print(length(which(apply(datSpecies,2,function (x) (effgroupe[i]-length(which(x[clusters==i]==0)))/effgroupe[i]*100)>90)))
#  print(length(which(apply(datSpecies,2,function (x) (clusters$size[i]-length(which(x[clusters$cluster==i]==0)))/clusters$size[i]*100)>50)))
#  print(length(which(apply(datSpecies,2,function (x) (clusters$size[i]-length(which(x[clusters$cluster==i]==0)))/clusters$size[i]*100)>90)))
  print(length(which(apply(cahacpclara$datSpecies,2,function (x) (cahacpclara$clusters$clusinfo[i,1]-length(which(x[cahacpclara$clusters$clustering==i]==0)))/cahacpclara$clusters$clusinfo[i,1]*100)>50)))
  print(length(which(apply(cahacpclara$datSpecies,2,function (x) (cahacpclara$clusters$clusinfo[i,1]-length(which(x[cahacpclara$clusters$clustering==i]==0)))/cahacpclara$clusters$clusinfo[i,1]*100)>90)))
}



# Test-values for the 5 bigger species (for test-value) (cluster by cluster)
# Valeurs-tests des 5 espèces qui ont les plus grandes valeurs-tests (classe par classe)
nbMainSpecies = min(5,dim(cahacpclara$target$tabnomespcib)[[2]])

namesVT=list()
for(i in 1:cahacpclara$nbClust){
  print("-----------------------------------------------------------------")
  print(paste("Cluster ",i))
  namesVT=rbind(namesVT,cahacpclara$target$tabnomespcib[i,1:nbMainSpecies])
  target5=cahacpclara$target$tabnomespcib[i,1:nbMainSpecies][!is.na(cahacpclara$target$tabnomespcib[i,1:nbMainSpecies])]
  print(target5)
  print(sapply(target5,function(x) cahacpclara$mprofil[i,x]))
  print(sapply(target5,function(x) cahacpclara$resval[x,i]))
}

# Percentage of mean catch for the 5 main caught species (cluster by cluster)
# Pourcentage de capture moyenne des 5 espèces les plus capturées (classe par classe)
nbMainSpecies = min(5,dim(cahacpclara$mprofil)[[2]])

dim(cahacpclara$mprofil)
namesCapt=list()
for(i in 1:cahacpclara$nbClust){
  print("-----------------------------------------------------------------")
  print(paste("Cluster ",i))
  namesCapt=rbind(namesCapt,names(cahacpclara$mprofil[i,][order(cahacpclara$mprofil[i,],decreasing=T)][1:nbMainSpecies]))
  print(cahacpclara$mprofil[i,][order(cahacpclara$mprofil[i,],decreasing=T)][1:nbMainSpecies])
  print(cahacpclara$resval[names(cahacpclara$mprofil[i,][order(cahacpclara$mprofil[i,],decreasing=T)][1:nbMainSpecies]),i])
}

# Percentage of cluster's logevents which catch these species (cluster by cluster)
# Pourcentage des logevents de la classe qui capturent ces espèces (classe par classe)
namesClust=cbind(namesVT,namesCapt)
for(i in 1:cahacpclara$nbClust){
  print("-----------------------------------------------------------------")
  print(paste("Cluster ",i))
  print(unique(unlist(namesClust[i,])))
#  print(unlist(lapply(as.list(unique(unlist(namesClust[i,]))),function(x) (sizeClusters[i]-length(which(datSpecies[clusters==i,x]==0)))/sizeClusters[i]*100)))
#  print(unlist(lapply(as.list(unique(unlist(namesClust[i,]))),function(x) (effgroupe[i]-length(which(datSpecies[clusters==i,x]==0)))/effgroupe[i]*100)))
#  print(unlist(lapply(as.list(unique(unlist(namesClust[i,]))),function(x) (clusters$size[i]-length(which(datSpecies[clusters$cluster==i,x]==0)))/clusters$size[i]*100)))
  print(unlist(lapply(as.list(unique(unlist(namesClust[i,]))),
  function(x) (cahacpclara$clusters$clusinfo[i,1]-length(which(cahacpclara$datSpecies[cahacpclara$clusters$clustering==i,x]==0)))/cahacpclara$clusters$clusinfo[i,1]*100)))
}

}
