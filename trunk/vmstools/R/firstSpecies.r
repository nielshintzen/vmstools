# SELECTION OF MAIN SPECIES THANKS TO THE ORDINATION METHOD

firstSpecies=function(dat, datLog, clusters){

  # Modification of the name of species
#  names(dat)[-1] <- unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))
#  names(dat)
  p=ncol(dat)
  firstSp=list()
  nbLogFirstSp=list()

  # Determine the first species for each logevent (the catchest species of the logevent)
  firstSp=apply(dat[,2:p],1,function(x) names(which.max(x)))
  
  # Calculate the number of logevent by first species (number of logevents by cluster)
  nbLogFirstSp=sapply(unique(unlist(firstSp)),function(x) table(firstSp==x)[2])

  names(nbLogFirstSp)=substr(names(nbLogFirstSp),1,3)
  
  # Nombre de logevents associés à chaque espèce principale issue de la méthode d'ordination (avec la 1ère espèce)
  png("Number of logevents by main species.png", width = 1200, height = 800)
  barplot(nbLogFirstSp[order(nbLogFirstSp,decreasing=T)], main="Number of logevents by main species", xlab="Main species", ylab="Number of Logevents", las=3)
  dev.off()
  # Nombre de logevents associés à chaque espèce principale issue de la méthode d'ordination (sans la 1ère espèce)
  png(paste("Number of logevents by main species without ",names(nbLogFirstSp[order(nbLogFirstSp,decreasing=T)])[1],".png",sep=""), width = 1200, height = 800)
  barplot(nbLogFirstSp[order(nbLogFirstSp,decreasing=T)][-1], main=paste("Number of logevents by main species without",names(nbLogFirstSp[order(nbLogFirstSp,decreasing=T)])[1],sep=" "), xlab="Main species", ylab="Number of Logevents", las=3)
  dev.off()
                                  
  
  # Nom des espèces principales issues de la méthode d'ordination
  nomEspOrdi=names(nbLogFirstSp)
  nbEspOrdi=length(nomEspOrdi)
  nomEspOrdi=nomEspOrdi[order(nomEspOrdi,decreasing=F)]
  # Initialisation
  numEspOrdiClarai=numeric()
  nomEspOrdiPlot=character()
  nbLogEspOrdiPlot=numeric()
  nbClust=length(unique(clusters))
  
  for(i in 1:nbClust){
    # nom des espèces principales issues d'ordination associées aux logevents de la classe CLARA i
    nomEspOrdiClarai=unique(unlist(firstSp[which(clusters==i)]))
    # leur numéro
    numEspOrdiClarai=sapply(nomEspOrdiClarai, function(x) which(nomEspOrdi==x))
    
    # nom des espèces principales à afficher sur le barplot de la classe i
    nomEspOrdiPloti=rep("",nbEspOrdi)
    nomEspOrdiPloti[numEspOrdiClarai]=nomEspOrdiClarai
    nomEspOrdiPlot=rbind(nomEspOrdiPlot,nomEspOrdiPloti)
    
    # pour la classe CLARA i, nombre de logevents de la classe CLARA i associé à chacune des espèces principales issues d'ordination
    nbLogEspOrdiPloti=rep(0,nbEspOrdi)
    nbLogEspOrdiPloti[numEspOrdiClarai]=sapply(nomEspOrdiClarai,function(x) rev(table(firstSp[which(clusters==i)]==x))[1])
    nbLogEspOrdiPlot=rbind(nbLogEspOrdiPlot,nbLogEspOrdiPloti)
  }
  colnames(nbLogEspOrdiPlot)=names(nbLogFirstSp)[order(names(nbLogFirstSp))]
  rownames(nbLogEspOrdiPlot)=paste("nbLogFirstSpCluster",seq=1:nbClust,sep="")
  
  
  # Comparaison 
  compNbLogClustVsFirstSp=table(clusters,firstSp)
  rownames(compNbLogClustVsFirstSp)=paste("Cluster",seq(1:nbClust),sep="")
  
  
#  pourcentNbLogEspOrdiPlot=t(apply(nbLogEspOrdiPlot,1,function(x) x/sum(x)*100))
#  barplot(pourcentNbLogEspOrdiPlot)
#  
#  write.table(pourcentNbLogEspOrdiPlot, file="pourcentNbLogEspOrdiPlot.txt", quote=T, dec='.', sep=';', col.names=T, row.names=T)
 
  
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
  plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
  plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
  plot(datLog[,1], datLog[,3], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
  plot(datLog[,1], datLog[,4], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
  plot(datLog[,2], datLog[,4], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
  plot(datLog[,3], datLog[,4], pch=21, bg=rainbow(length(nbLogFirstSp))[as.numeric(numFirstSp)], main="Projection of the classification by ordination on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
  par(op)
  dev.off()

  return(list(firstSp=firstSp, nbLogFirstSp=nbLogFirstSp, compNbLogClustVsFirstSp=compNbLogClustVsFirstSp))
}