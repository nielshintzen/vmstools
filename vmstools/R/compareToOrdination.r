###########################################################################
#   DETERMINE MÉTIERS LEVEL 7 AND 5 WITH THE ORDINATION METHOD            #
#   (FIRST SPECIES IN CATCH OR FIRST GROUP OF SPECIES IN CATCH)           #
#                                                                         #
#   DETERMINE MÉTIERS LEVEL5 FROM MÉTIERS LEVEL7 FOUND WITH THE           #
#   MULTIVARIATE CLASSIFICATION (DETERMINING THE CHARACTERISTIC           #
#   SPECIES FOR EACH CLUSTER)                                             #
#                                                                         #
#   COMPARE THE MÉTIERS FOUND WITH ORDINATION METHOD AND                  #
#   WITH MULTIVARIATE CLASSIFICATION                                      #
###########################################################################  
  
compareToOrdination=function(dat, Step2, clusters){

  ##########################
  # FIRST SPECIES IN CATCH #
  ##########################
  
  # Modification of the name of species
  #  names(dat)[-1] <- unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))
  #  names(dat)
  dat=as.matrix(dat[,2:ncol(dat)])
  p=ncol(dat)
  firstSp=list()
  nbLogFirstSp=list()

  # Determine the first species for each logevent (the catchest species of the logevent)
  firstSp=apply(dat,1,function(x) names(which.max(x)))
  
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
  #firstSp=sapply(nomEspOrdi, function(x) firstSp[which(firstSp==x)]=espOrdi[2,which(espOrdi[1,]==x)])
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
  
  # Comparison metiers L7 (clusters clara_step3) vs metiers L7 (first species in catch)
  compNbLogClustVsFirstSp=table(clusters,firstSp)
  rownames(compNbLogClustVsFirstSp)=paste("Cluster",seq(1:nbClust),sep="")
  
  # Determine the metier in level 5 from the metier in level 7 based on the first species in catch (logevent by logevent)
  logSpeciesL5=unlist(sapply(firstSp,level5))      ## pb espèces sans correspondance vers groupe d'espèces
#  unique(logSpeciesL5)
#  barplot(table(logSpeciesL5), main="Number of logevents by group of species level 5", xlab="Groups of species level 5", ylab="Number of logevents")
 
 
   
  
  ###################################
  # FIRST GROUP OF SPECIES IN CATCH #
  ###################################
  
  # Determine the metier of each logevent thanks to the first group of species in catch (level 5) of the logevent
  datL5=dat
  colnames(datL5)=lapply(as.list(colnames(datL5)), level7to5)   ## quand on aura toutes les correspondances espèces vers groupes d'espèces, on pourra utiliser level5 à la place de level7to5.
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
  logGroupL5=unique(colnames(datL5))[apply(tab,1,which.max)]
  
  
  # Compare the metiers (level 5) found thanks to the first species in catch vs the first group of species in catch
  #compL5SpeciesVsGroup=table(logSpeciesL5,logGroupL5)
  
  
  
  ###############################################################
  # LEVELS 7 AND 5 FOR METIERS FROM MULTIVARIATE CLASSIFICATION #
  ###############################################################
  
  # We are taking species of the cluster (by decreasing order of catch) until having at least seuilCatch% of the cluster's total catch
  seuilCatch=75
  nameSpecies=colnames(dat)
  
  # Total quantity caught species by species by cluster
  sumColByClust=matrix(0,nrow=nbClust,ncol=p)
  for(cl in 1:nbClust){
    for(i in 1:p) sumColByClust[cl,i]=sum(dat[which(clusters==cl),i])
  }
  colnames(sumColByClust)=nameSpecies
  
  
  # Percent of each species in the cluster's total catch
  propEspByCluster=matrix(0,nrow=nbClust,ncol=p)
  colnames(propEspByCluster)=nameSpecies
  for(cl in 1:nbClust) propEspByCluster[cl,]=sumColByClust[cl,]/sum(sumColByClust[cl,])*100
  
  # Columns number of each species by decreasing order of catch
  numEspByCluster=matrix(0,nrow=nbClust,ncol=p)
  for(cl in 1:nbClust) numEspByCluster[cl,]=order(propEspByCluster[cl,],decreasing=T)
  # Cumulative percentage of species catch in the cluster's total catch
  propEspByClusterCum=matrix(0,nrow=nbClust,ncol=p)
  for(cl in 1:nbClust)  propEspByClusterCum[cl,]=cumsum(sort(propEspByCluster[cl,],decreasing=T))
  
   
  espSelByCluster=matrix(0,nrow=nbClust,ncol=p)
  nomEspSelByCluster=matrix(NA,nrow=nbClust,ncol=p)
  for(cl in 1:nbClust){ 
    pourcent=which(propEspByClusterCum[cl,]<=seuilCatch)
    # Name of selected species by cluster
    espSelByCluster[cl,1:(length(pourcent)+1)]=numEspByCluster[cl,1:(length(pourcent)+1)]
    nomEspSelByCluster[cl,1:(length(pourcent)+1)]=nameSpecies[espSelByCluster[cl,1:(length(pourcent)+1)]]
  }
      
  # Maximum number of species in the table
  maxColNomEspSelByCluster=max(unlist(lapply((apply(!is.na(nomEspSelByCluster),1,which)),length)))
  
  nomEspSelByCluster=nomEspSelByCluster[,1:maxColNomEspSelByCluster]
  # List of selected species names by cluster (= List of metiers level 7)
  listNomEspSelByCluster=list()
  for(cl in 1:nbClust) listNomEspSelByCluster[[cl]]=unlist(nomEspSelByCluster[cl,which(!is.na(nomEspSelByCluster[cl,]))])
  
  # List of metiers level 5.7
  listMetiersLevel57=lapply(listNomEspSelByCluster, level5)
  
  # Metiers (level 5.7) of each logevent found thanks to the getMetierClusters method.
  clustersStep3L5=unlist(sapply(clusters,function(x) listMetiersLevel57[x]))
  
  
  #################
  #  COMPARISONS  # 
  #################
  
  # Compare the metiers of getMetierClusters (level 5.7) vs the first group of species in catch (level 5)
  #compL5ClustStep3VsGroup=table(clustersStep3L5,logGroupL5)  !! pb regroupement des clusters qui correspondent au même groupe d'espèce niveau 5.
  clustStep3L5=clusters
  compL57vsL5=table(clustStep3L5,logGroupL5)
  rownames(compL57vsL5)=unlist(listMetiersLevel57)
  
  # Compare the metiers (level 5.7) vs metiers (level 7)
  clustStep3L7=clusters
  compL57vsL7=table(clustersStep3L5,clustStep3L7)
  colnames(compL57vsL7)=paste("Clust",seq=1:nbClust)
  
  # Compare the metiers (level 5) vs metiers (level 7)
  compL5vsL7=table(logGroupL5,clustStep3L7)
  colnames(compL5vsL7)=paste("Clust",seq=1:nbClust)
  
  # Compare the metiers (level 5.7) found thanks to the getMetierClusters vs the first species in catch (level 5)
  #compL5ClustStep3VsSpecies=table(clustersStep3L5,logSpeciesL5)    !! pb regroupement des clusters qui correspondent au même groupe d'espèce niveau 5.
#  compL5ClustStep3VsSpecies=table(clustStep3L5,logSpeciesL5)         ## pb avec logSpeciesL5
#  rownames(compL5ClustStep3VsSpecies)=unlist(listMetiersLevel57)
  

  
  
  
  
#  pourcentNbLogEspOrdiPlot=t(apply(nbLogEspOrdiPlot,1,function(x) x/sum(x)*100))
#  barplot(pourcentNbLogEspOrdiPlot)
#  
#  write.table(pourcentNbLogEspOrdiPlot, file="pourcentNbLogEspOrdiPlot.txt", quote=T, dec='.', sep=';', col.names=T, row.names=T)

  tabCompL57vsL5=as.data.frame(apply(format(compL57vsL5), 2, as.numeric))
  tabCompL57vsL7=as.data.frame(apply(format(compL57vsL7), 2, as.numeric))
  tabCompL5vsL7=as.data.frame(apply(format(compL5vsL7), 2, as.numeric))

  db="compClassifvsOrdin.xls"
  if(db %in% list.files()) file.remove(db)
  channel <- odbcConnectExcel(xls.file = db,readOnly=FALSE)
  sqlSave(channel, cbind(rownames(compL57vsL5),tabCompL57vsL5), tablename = "compL57vsL5")
  sqlSave(channel, cbind(rownames(compL57vsL7),tabCompL57vsL7), tablename = "compL57vsL7")
  sqlSave(channel, cbind(rownames(compL5vsL7),tabCompL5vsL7), tablename = "compL5vsL7")
  odbcClose(channel)
  
  

  return(list(nbLogFirstSp=nbLogFirstSp, compNbLogClustVsFirstSp=compNbLogClustVsFirstSp,
   compL57vsL5=compL57vsL5, compL57vsL7=compL57vsL7, compL5vsL7=compL5vsL7))
}