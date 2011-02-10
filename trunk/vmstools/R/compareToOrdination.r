###########################################################################
#   DETERMINE MÉTIERS LEVEL 7 AND 5 WITH TWO ORDINATION METHODS           #
#   (FIRST SPECIES IN CATCH AND FIRST GROUP OF SPECIES IN CATCH)          #
#                                                                         #
#   DETERMINE MÉTIERS LEVEL5 FROM MÉTIERS LEVEL7 FOUND WITH THE           #
#   MULTIVARIATE CLASSIFICATION (DETERMINING THE TARGET SPECIES           #
#   FOR EACH CLUSTER)                                                     #
#                                                                         #
#   COMPARE THE MÉTIERS FOUND WITH THE ORDINATION METHODS AND             #
#   WITH THE MULTIVARIATE CLASSIFICATION                                  #
###########################################################################  
  
compareToOrdination=function(dat, Step2, clusters, tabClusters){

  # Load the table linking 3A-CODE (FAO CODE of species) to the species assemblage (level 5).
  data(correspLevel7to5)
  
  ##########################
  # FIRST SPECIES IN CATCH #
  ##########################
  
  dat=as.matrix(dat[,2:ncol(dat)])
  p=ncol(dat)
  firstSp=list()
  nbLogFirstSp=list()

  # Determine the first species for each logevent (the catchest species of the logevent)
  firstSp=apply(dat,1,function(x) names(which.max(x)))
  metiersFirstSpeciesL7=firstSp
  
  # Calculate the number of logevent by first species (number of logevents by cluster)
  nbLogFirstSp=sapply(unique(unlist(firstSp)),function(x) table(firstSp==x)[2])

  names(nbLogFirstSp)=substr(names(nbLogFirstSp),1,3)
  
  # Number of logevents linked to each species from the method "first species in catch"
  png("Number of logevents by main species.png", width = 1200, height = 800)
  barplot(nbLogFirstSp[order(nbLogFirstSp,decreasing=T)], main="Number of logevents by main species", xlab="Main species", ylab="Number of Logevents", las=3)
  dev.off()
  # idem without the species with the biggest number of logevents
  png(paste("Number of logevents by main species without ",names(nbLogFirstSp[order(nbLogFirstSp,decreasing=T)])[1],".png",sep=""), width = 1200, height = 800)
  barplot(nbLogFirstSp[order(nbLogFirstSp,decreasing=T)][-1], main=paste("Number of logevents by main species without",names(nbLogFirstSp[order(nbLogFirstSp,decreasing=T)])[1],sep=" "), xlab="Main species", ylab="Number of Logevents", las=3)
  dev.off()
                                  
  
  # Names of species from the ordination method
  nomEspOrdi=names(nbLogFirstSp)
  nbEspOrdi=length(nomEspOrdi)
  nomEspOrdi=nomEspOrdi[order(nomEspOrdi,decreasing=F)]
  # Initialization
  numEspOrdiClarai=numeric()
  nomEspOrdiPlot=character()
  nbLogEspOrdiPlot=numeric()
  nbClust=length(unique(clusters))
  
  # For each cluster Clara
  for(i in 1:nbClust){
    # Names of species from ordination method linked to the logevents of cluster Clara i
    nomEspOrdiClarai=unique(unlist(firstSp[which(clusters==i)]))
    # their number
    numEspOrdiClarai=sapply(nomEspOrdiClarai, function(x) which(nomEspOrdi==x))
    
    # Names of species to plot in the barplot of cluster i
    nomEspOrdiPloti=rep("",nbEspOrdi)
    nomEspOrdiPloti[numEspOrdiClarai]=nomEspOrdiClarai
    nomEspOrdiPlot=rbind(nomEspOrdiPlot,nomEspOrdiPloti)
    
    # for the cluster i, number of logevents of cluster i linked to each species from ordination method
    nbLogEspOrdiPloti=rep(0,nbEspOrdi)
    nbLogEspOrdiPloti[numEspOrdiClarai]=sapply(nomEspOrdiClarai,function(x) rev(table(firstSp[which(clusters==i)]==x))[1])
    nbLogEspOrdiPlot=rbind(nbLogEspOrdiPlot,nbLogEspOrdiPloti)
  }
  colnames(nbLogEspOrdiPlot)=names(nbLogFirstSp)[order(names(nbLogFirstSp))]
  rownames(nbLogEspOrdiPlot)=paste("nbLogFirstSpCluster",seq=1:nbClust,sep="")
  
  # Graphics
  png("Number of logevents linked to each ordination's main species by CLARA's cluster.png", width = 1200, height = 800)
  op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
  for(i in 1:nbClust){  
    op2 <- par(las=2)
    barplot(nbLogEspOrdiPlot[i,],names.arg=nomEspOrdiPlot[i,], las=3)
    par(op2)
    mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
  }
  par(op)
  title(main=paste("Number of logevents linked to each ordination's main species by CLARA's cluster","\n","\n",sep="")) 
  dev.off()
  
  # Projections on the first factorial plans
  nomEspOrdi=substr(names(nbLogFirstSp),1,3)
  espOrdi=matrix()
  espOrdi=rbind(names(nbLogFirstSp),seq(1:length(nbLogFirstSp)))
  numFirstSp=firstSp
  for(i in 1:length(nbLogFirstSp)){
    numFirstSp[which(firstSp==nomEspOrdi[i])]=espOrdi[2,which(espOrdi[1,]==nomEspOrdi[i])]
  }

  # Projections on the first factorial plans
  png("Ordination_Projections.png", width = 1200, height = 800)
  op <- par(mfrow=c(2,3))
  plot(Step2[,1], Step2[,2], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
  plot(Step2[,2], Step2[,3], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
  plot(Step2[,1], Step2[,3], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
  plot(Step2[,1], Step2[,4], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
  plot(Step2[,2], Step2[,4], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
  plot(Step2[,3], Step2[,4], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
  par(op)
  dev.off()
  
  
  # Determine the metier in level 5 from the metier in level 7 based on the first species in catch (logevent by logevent)
#  metiersFirstSpeciesL5=unlist(sapply(firstSp,level5))      ## pb espèces sans correspondance vers groupe d'espèces
#  unique(metiersFirstSpeciesL5)
#  barplot(table(metiersFirstSpeciesL5), main="Number of logevents by group of species level 5", xlab="Groups of species level 5", ylab="Number of logevents")
 
  print("Please, be patient...") 
  metiersFirstSpeciesL5=lapply(as.list(firstSp), function(x) if(length(which(correspLevel7to5[,"X3A_CODE"]==x))==0){print(paste(x," : unknown species, classed in \'FIF\' group",sep=""));"FIF"} 
                                                    else correspLevel7to5[which(correspLevel7to5[,"X3A_CODE"]==x),"DCF_species_level5_COD"])
  metiersFirstSpeciesL5=unlist(lapply(metiersFirstSpeciesL5, function(x) as.character(x)))
   
   
   
                                      
  ###################################
  # FIRST GROUP OF SPECIES IN CATCH #
  ###################################
  
  # Determine the metier of each logevent thanks to the first group of species in catch (level 5) of the logevent
  datL5=dat
  #colnames(datL5)=lapply(as.list(colnames(datL5)), level7to5)   ## quand on aura toutes les correspondances espèces vers groupes d'espèces, on pourra utiliser level5 à la place de level7to5.
  groupColSpecies=lapply(as.list(colnames(datL5)), function(x) if(length(which(correspLevel7to5[,"X3A_CODE"]==x))==0){print(paste(x," : unknown species, classed in \'FIF\' group",sep=""));"FIF"} 
                                                               else correspLevel7to5[which(correspLevel7to5[,"X3A_CODE"]==x),"DCF_species_level5_COD"])  
  groupColSpecies=unlist(lapply(groupColSpecies, function(x) as.character(x)))
  colnames(datL5)=groupColSpecies
  nbC5=length(unique(colnames(datL5)))
  tab=matrix(NA,nrow=nrow(datL5),ncol=nbC5)
  colnames(tab)=unique(colnames(datL5))
  for(i in 1:nbC5){
    # Sub-data containing only one metier level 5.
    subdati=datL5[,which(colnames(datL5)==(unique(colnames(datL5))[i]))]
    if(length(which(colnames(datL5)==(unique(colnames(datL5))[i])))==1) tab[,i]=subdati
    else
      tab[,i]=apply(subdati,1,sum)
  }
  # Metiers (level 5) of each logevent found thanks to the first group of species in catch (level 5)
  metiersFirstGroupL5=unique(colnames(datL5))[apply(tab,1,which.max)]
  
  


  
  
  ###############################################################
  # LEVELS 7 AND 5 FOR METIERS FROM MULTIVARIATE CLASSIFICATION #
  ###############################################################
  
  # We are taking all species in the tabClusters until having at least seuilCatch% of cumulated "% Catch".
  seuilCatch=75
  
  targetSpeciesByCluster=matrix(NA,nrow=nbClust,ncol=10)
  for(i in 1:nbClust){
    percCatchCum=cumsum(as.numeric(tabClusters[,"% Catch",i]))
    nbSpSel=length(which(percCatchCum<seuilCatch))+1
    targetSpeciesByCluster[i,1:nbSpSel]=tabClusters[1:nbSpSel,"FAO",i]
  }
  
  # Maximum number of species in the table
  maxColNomEspSelByCluster=max(unlist(lapply((apply(!is.na(targetSpeciesByCluster),1,which)),length)))
  targetSpeciesByCluster=targetSpeciesByCluster[,1:maxColNomEspSelByCluster]
  
  # List of selected species names by cluster (= List of metiers level 7)
  listTargetSpeciesByCluster=list()
  for(cl in 1:nbClust) listTargetSpeciesByCluster[[cl]]=unlist(targetSpeciesByCluster[cl,which(!is.na(targetSpeciesByCluster[cl,]))])
  
  # List of metiers level 57
  #listMetiersLevel57=lapply(listTargetSpeciesByCluster, level5)
  listMetiersLevel57=list()
  for (i in 1:nbClust){
    metiersClusteri=lapply(listTargetSpeciesByCluster[[i]], function(x) if(length(which(correspLevel7to5[,"X3A_CODE"]==x))==0){print(paste(x," : unknown species, classed in \'FIF\' group",sep=""));"FIF"}
                                                        else correspLevel7to5[which(correspLevel7to5[,"X3A_CODE"]==x),"DCF_species_level5_COD"])
    metiersClusteri=as.character(unique(unlist(metiersClusteri)))
    metiersClusteri=paste(unlist(strsplit(metiersClusteri," ")),collapse=" ")
    listMetiersLevel57[[i]]=metiersClusteri
  }
  
  
  # Metiers (level 57) of each logevent found thanks to the getMetierClusters method.
  metiersClustersL5=unlist(sapply(clusters,function(x) listMetiersLevel57[x]))
  

  
  
  
     
   ## Changement dans la définition des espèces cibles (niveau 7) : on ne prend plus les espèces les plus capturées dans la classe
   ## jusqu'à avoir au moins 75% de la capture totale de la classe
   ## on prend maintenant les espèces présentes dans le tableau de synthèse de la classe jusqu'à avoir un %Catch cumulé au moins 
   ## égal à 75% 
   ## (Rappel : %Catch d'une espèce = la moyenne des pourcentages de capture de cette espèce pour tous les logevents de la classe)
  
#  nameSpecies=colnames(dat)
#  
#  # Total quantity caught species by species by cluster
#  sumColByClust=matrix(0,nrow=nbClust,ncol=p)
#  for(cl in 1:nbClust){
#    for(i in 1:p) sumColByClust[cl,i]=sum(dat[which(clusters==cl),i])
#  }
#  colnames(sumColByClust)=nameSpecies
#  
#  
#  # Percent of each species in the cluster's total catch
#  propEspByCluster=matrix(0,nrow=nbClust,ncol=p)
#  colnames(propEspByCluster)=nameSpecies
#  for(cl in 1:nbClust) propEspByCluster[cl,]=sumColByClust[cl,]/sum(sumColByClust[cl,])*100
#  
#  # Columns number of each species by decreasing order of catch
#  numEspByCluster=matrix(0,nrow=nbClust,ncol=p)
#  for(cl in 1:nbClust) numEspByCluster[cl,]=order(propEspByCluster[cl,],decreasing=T)
#  # Cumulative percentage of species catch in the cluster's total catch
#  propEspByClusterCum=matrix(0,nrow=nbClust,ncol=p)
#  for(cl in 1:nbClust)  propEspByClusterCum[cl,]=cumsum(sort(propEspByCluster[cl,],decreasing=T))
#  
#   
#  espSelByCluster=matrix(0,nrow=nbClust,ncol=p)
#  nomEspSelByCluster=matrix(NA,nrow=nbClust,ncol=p)
#  for(cl in 1:nbClust){ 
#    pourcent=which(propEspByClusterCum[cl,]<seuilCatch)
#    # Name of selected species by cluster
#    espSelByCluster[cl,1:(length(pourcent)+1)]=numEspByCluster[cl,1:(length(pourcent)+1)]
#    nomEspSelByCluster[cl,1:(length(pourcent)+1)]=nameSpecies[espSelByCluster[cl,1:(length(pourcent)+1)]]
#  }
#      
#  # Maximum number of species in the table
#  maxColNomEspSelByCluster=max(unlist(lapply((apply(!is.na(nomEspSelByCluster),1,which)),length)))
#  
#  nomEspSelByCluster=nomEspSelByCluster[,1:maxColNomEspSelByCluster]
#  # List of selected species names by cluster (= List of metiers level 7)
#  listNomEspSelByCluster=list()
#  for(cl in 1:nbClust) listNomEspSelByCluster[[cl]]=unlist(nomEspSelByCluster[cl,which(!is.na(nomEspSelByCluster[cl,]))])
#  
#  # List of metiers level 5.7
#  listMetiersLevel57=lapply(listNomEspSelByCluster, level5)
#  
#  # Metiers (level 5.7) of each logevent found thanks to the getMetierClusters method.
#  clustersStep3L5=unlist(sapply(clusters,function(x) listMetiersLevel57[x]))
  
  
  
  #################
  #  COMPARISONS  # 
  #################
  
  # Compare the metiers Clusters (from getMetierClusters) (level 7) vs the metier First Species in catch (level 7)
  compClustersL7vsFirstSpeciesL7=table(clusters,metiersFirstSpeciesL7)
  rownames(compClustersL7vsFirstSpeciesL7)=paste("Clust",seq=1:nbClust,sep="")
  
  # Compare the metiers Clusters (level 5.7) vs metiers First Species (level 5)
  clustStep3L5=clusters
  compClustersL5vsFirstSpeciesL5=table(clustStep3L5,metiersFirstSpeciesL5)    
  rownames(compClustersL5vsFirstSpeciesL5)=unlist(listMetiersLevel57) 
  
  # Compare the metiers Clusters (level 5.7) vs metiers First Group (level 5)
  #compL5ClustStep3VsGroup=table(clustersStep3L5,metiersFirstGroupL5)  !! pb regroupement des clusters qui correspondent au même groupe d'espèce niveau 5.
  clustStep3L5=clusters
  compClustersL5vsFirstGroupL5=table(clustStep3L5,metiersFirstGroupL5)
  rownames(compClustersL5vsFirstGroupL5)=unlist(listMetiersLevel57)
  
  # Compare the metiers First Species (level 5) vs First Group (level 5)
  compFirstSpeciesL5vsFirstGroupL5=table(metiersFirstSpeciesL5,metiersFirstGroupL5)
  
  
  # Less attractive
  # Compare the metiers Clusters (level 5.7) vs metiers Clusters (level 7)
  clustStep3L7=clusters
  compClustersL5vsClustersL7=table(metiersClustersL5,clustStep3L7)
  colnames(compClustersL5vsClustersL7)=paste("Clust",seq=1:nbClust)
  
  
           
  # Create csv tables
  dfTables=data.frame()
  dfTables=c("compClustersL7vsFirstSpeciesL7","compClustersL5vsFirstSpeciesL5",
  "compClustersL5vsFirstGroupL5","compFirstSpeciesL5vsFirstGroupL5","compClustersL5vsClustersL7")
  write.table(dfTables[1],file="tablesCompToOrdination.csv",append=TRUE,col.names=NA)
  write.table(compClustersL7vsFirstSpeciesL7,file="tablesCompToOrdination.csv",append=TRUE,col.names=NA)
  write.table(dfTables[2],file="tablesCompToOrdination.csv",append=TRUE,col.names=NA)
  write.table(compClustersL5vsFirstSpeciesL5,file="tablesCompToOrdination.csv",append=TRUE,col.names=NA)
  write.table(dfTables[3],file="tablesCompToOrdination.csv",append=TRUE,col.names=NA)
  write.table(compClustersL5vsFirstGroupL5,file="tablesCompToOrdination.csv",append=TRUE,col.names=NA)
  write.table(dfTables[4],file="tablesCompToOrdination.csv",append=TRUE,col.names=NA)
  write.table(compFirstSpeciesL5vsFirstGroupL5,file="tablesCompToOrdination.csv",append=TRUE,col.names=NA)
  write.table(dfTables[5],file="tablesCompToOrdination.csv",append=TRUE,col.names=NA)
  write.table(compClustersL5vsClustersL7,file="tablesCompToOrdination.csv",append=TRUE,col.names=NA)
      
  print("Done.")
  
  return(list(nbLogFirstSp=nbLogFirstSp, 
    compClustersL7vsFirstSpeciesL7=compClustersL7vsFirstSpeciesL7, 
    compClustersL5vsFirstSpeciesL5=compClustersL5vsFirstSpeciesL5,
    compClustersL5vsFirstGroupL5=compClustersL5vsFirstGroupL5,
    compFirstSpeciesL5vsFirstGroupL5=compFirstSpeciesL5vsFirstGroupL5,
    compClustersL5vsClustersL7=compClustersL5vsClustersL7))
}