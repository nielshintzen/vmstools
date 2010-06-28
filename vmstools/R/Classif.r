##############################################################################################################
#                                                                                                            #
#     CLASSIF : FUNCTION FOR SPECIES SELECTION, PCA, CLASSIFICATION OF LOGEVENTS AND MÉTIERS COMPUTATION     #
#                                                                                                            #
##############################################################################################################

                                                                                
#############################################################
#  PARAMETERS TO FILL IN DEPENDING OF THE SELECTED METHODS  #
#############################################################
#
# STEP 1 : Choose the method for the selection of species : 3 methods : CAH, total, logevent
#
# CAH : we view as main species, all species which are not belonging to the cluster with the smaller mean
# 2 parameters to fill in : param1 for the distance between individuals ("euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski")
#                           param2 for the agregative distance (distance between clusters) ("ward", "single", "complete", "average", "mcquitty", "median" or "centroid")
#
# total : we take all species by decreasing importance order until having a given percentage (70, 80, 90%) of the total capture
# 1 parameter to fill in : param1 = the threshold (70, 80, 90...)
#
# logevent : we take all species which represents at least a given percentage (10, 20, 30%) for at least a logevent
# 1 parameter to fill in : param1 = the threshold (10, 20, 30...)
#
#
# STEP 2 : PCA or not ?
#
# pca : we execute a PCA on the individuals "logevents", and compute the classification on the factorial coordinates of the logevents
# nopca : we choose to work with the raw data to compute the classification of logevents
#
#
# STEP 3 : Choose the method for the classification of logevents : CAH, KMEANS, PAM, CLARA ?
#
# CAH, 2 parameters to fill in : param3 for the distance between individuals ("euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski")
#                                param4 for the agregative distance (distance between clusters) ("ward", "single", "complete", "average", "mcquitty", "median" ou "centroid")
#
# KMEANS, no parameter to fill in
#
# PAM, no parameter to fill in
#
# CLARA, 1 parameter to fill in : param3 for the dissimilarity between individuals ("euclidean" or "manhattan")




# Classif : function for the selection of species (hac, total, logevent), PCA (pca, nopca), classification of logevents (hac, kmeans, pam, and clara) and métiers computation (thanks to test-values)
classif=function(dat,methSpecies="hac",param1="euclidean",param2="ward",pcaYesNo="pca",methMetier="clara",param3="euclidean",param4=NULL){

#first siplify the names of columns 
names(dat)[-1] <- unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))
names(dat)

# STEP 1 : SELECTION OF SPECIES

print("######## STEP 1 SELECTION OF MAIN SPECIES ########")

  if(methSpecies=="hac"){

    toutfait=FALSE
    p=length(dat[1,])   # Number of species
    # Transform quantities to proportions of total quantity caught by logevent
    propdat=transformation_proportion(dat[,2:p])
    names(propdat)
    # Transposing data
    table_var=table_variables(propdat)

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   

    # CAH
    esp_dist=dist(table_var, method=param1)
    cah_var=hclust(esp_dist, method=param2)

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   

    # Select the number of clusters by scree-test
    inerties.vector=cah_var$height[order(cah_var$height,decreasing=T)]
    nb.finalclusters=which(scree(inerties.vector)$epsilon<0)[1]

    # Dendogram cutting at the selected level
    cah_cluster_var=cutree(cah_var,k=nb.finalclusters)

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])  
                                                        
    # Selection of main species
    espprinc=select_species(dat[,2:p],cah_cluster_var)
    espdata=names(dat[,2:p])


    while(toutfait==FALSE){
      ind_princ=which(is.element(espdata,espprinc))
      ind_autres=setdiff(espdata,espprinc)
      princ=table_var[ind_princ,]
      resi=table_var[ind_autres,]

      # CAH on residual species group
      resiesp_dist=dist(resi, method=param1)
      resicah_var=hclust(resiesp_dist, method=param2)
      # Select the number of clusters by scree-test
      resiinerties.vector=resicah_var$height[order(resicah_var$height,decreasing=T)]
      resinbclusters=which(scree(resiinerties.vector)$epsilon<0)[1]
      # Cut the dendogram at the selected level
      resicah_cluster_var=cutree(resicah_var,k=resinbclusters)

      # Number and names of species by cluster, mean of capture by species for each cluster
      long=numeric()
      somclus=numeric()
      moyclus=numeric()
      for(i in 1:resinbclusters){
        long[i]=length(which(resicah_cluster_var==i))
        nomclus=names(which(resicah_cluster_var==i))
        somclus[i]=sum(resi[nomclus,])
        moyclus[i]=somclus[i]/long[i]
      }
      # we are taking the cluster which represents the greatest part of capture
      indclusmax=which.max(moyclus)
      names(which(resicah_cluster_var==indclusmax))
      # if the cluster owns only one species, we are adding the species to the main species and restart the CAH with residuals species
      if(length(which(resicah_cluster_var==indclusmax))==1){
        espprinc=c(espprinc,names(which(resicah_cluster_var==indclusmax)))
        toutfait=FALSE
      }else{   # else we are stopping here and keeping the main species of the start
        toutfait=TRUE
      }

    }

    # Aggregate species of the cluster with the smaller mean in variable AUTRES (i.e. OtherSpecies)
    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])

    datSpecies=building_tab_pca(propdat,espprinc)

  }

  if (length(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
  }              
  
  
  if(methSpecies=="total"){

    propdat=transformation_proportion(dat[,2:p])
    p=length(dat[1,])   # Number of species
    # Total quantity caught by species
    sumcol=numeric()
    for(i in 2:p){
      sumcol[i]=sum(dat[,i])
    }
    #clu 
    # sumcol <- apply(dat[,-1],2,sum)
    
    # Total quantity caught
    sumtotale=sum(sumcol,na.rm=T)
    # Percent of each species in the total catch
    propesp=sumcol/sumtotale*100
    # by decreasing order
    propespdec=propesp[order(propesp,decreasing=T)]
    # Columns number of each species by decreasing order of capture
    numesp=order(propesp,decreasing=T)
    # Cumulative percent of catch
    propespcum=cumsum(propespdec)
    # We are taking all species until having at least param1% of total catch
    seuil=param1
    pourcent=which(propespcum<seuil)
    espsel=numesp[1:(length(pourcent)+1)]
    # We are taking the name of selected species
    nomespsel=names(dat[espsel])
    # We are bulding the table with main species and aggregated other species
    datSpecies=building_tab_pca(propdat,nomespsel)

  }
  if (length(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
  }


  if(methSpecies=="logevent"){

    p=length(dat[1,])   # Number of species
    # Transform quantities to proportions of total quantity caught by logevent
    propdat=transformation_proportion(dat[,2:p])
 #   n=length(propdat[,1])       # Number of logevents
 #   nomespsel=numeric()
    #seuil=seuillogevent(propdat)
    seuil=param1
    # Selection of species making up over param1% of logevent's captures
    #  clu
    pourcent <- apply(propdat,1,function(x) which(x>seuil))   
    nomespsel <- names(propdat)[unique(unlist(pourcent))]  
#    for(i in 1:n){
#      # We are sorting the logevent by caught species decreasing order
#      logeventdec=propdat[i,][order(propdat[i,],decreasing=T)]
#      # We are recording the number of these species
#      numesp=order(propdat[i,],decreasing=T)
#      # We are selecting species making up over param1% of logevent's captures
#      pourcent=which(logeventdec>seuil)
#      # if some species are selected, we are recording their name
#      if(length(pourcent)>0){
#        espsel=numesp[1:length(pourcent)]
#        nomespsel[i]=names(dat[,2:p][espsel])
#      }else{      # else we are recording NA
#        espsel=NULL
#        nomespsel[i]="NA"
#      }
#    }
    # We are bulding the table with main species and aggregated other species
    #datSpecies=building_tab_pca(dat[,2:p],nomespsel)
    datSpecies=building_tab_pca(propdat,nomespsel)

  }
  if (length(objects()[-which(objects() %in% c('methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
  }



# STEP 2 : PCA OR NOT

print("######## STEP 2 PCA/NO PCA ON CATCH PROFILES ########")

  if(pcaYesNo=="pca"){
    # PCA (Principal Component Analysis)
    log.pca <- PCA(datSpecies, graph=T)
    
    savePlot(filename='Species projection on the two first factorial axis', type='png', restoreConsole = TRUE)
    dev.off()
    savePlot(filename='Individuals projection on the two first factorial axis', type='png', restoreConsole = TRUE)
    dev.off()
    
    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   

    graphics.off()

    # Eigenvalues and relative graphics
    log.pca$eig                                                  

    png("Eigen values.png", width = 1200, height = 800)
    x=1:length(log.pca$eig[,1])
    barplot(log.pca$eig[,1],names.arg=x, main="Eigen values")
    dev.off()
    png("Percentage of Inertia.png", width = 1200, height = 800)
    barplot(log.pca$eig[,2],names.arg=x, main="Percentage of Inertia of factorial axis", xlab="Axis", ylab="% of Inertia")
    dev.off()
    png("Cumulative Percentage of Inertia.png", width = 1200, height = 800) 
    color=rep("grey",length(log.pca$eig[,1]))
    numaxe=min(which(log.pca$eig[,3]>70))
    color[1:numaxe]="green"
    barplot(log.pca$eig[,3],names.arg=x, col=color, main="Cumulative Percentage of Inertia of factorial axis", xlab="Axis", ylab="% of Inertia")
    abline(h=70, col="red")
    text(1,72, "70% of Inertia", col = "red", adj = c(0, -.1))
    dev.off()   
    
    # Data frame given eigenvalues, inertia and cumulative inertia of factorial axis
    tabInertia=data.frame(cbind(Axis=1:length(log.pca$eig[,1]), Eigenvalues=log.pca$eig[,1], Inertia=log.pca$eig[,2], CumulativeInertia=log.pca$eig[,3]))                  

    # Projection of variables Species on the first factorial axis
    png("Projection of Species on first factorial axis.png", width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(log.pca,choix="var",axes = c(1, 2),new.plot=FALSE,lim.cos2.var = 0.5)
    plot(log.pca,choix="var",axes = c(2, 3),new.plot=FALSE,lim.cos2.var = 0.5)
    plot(log.pca,choix="var",axes = c(1, 3),new.plot=FALSE,lim.cos2.var = 0.5)
    plot(log.pca,choix="var",axes = c(1, 4),new.plot=FALSE,lim.cos2.var = 0.5)
    plot(log.pca,choix="var",axes = c(2, 4),new.plot=FALSE,lim.cos2.var = 0.5)
    plot(log.pca,choix="var",axes = c(3, 4),new.plot=FALSE,lim.cos2.var = 0.5)
    par(op)
    dev.off() 


    # Determine the number of axis
    nbaxes=which(log.pca$eig[,3]>70)[1]   # we are taking the axis until having 70% of total inertia
    log.coa=PCA(datSpecies, graph=F, ncp=nbaxes)
    options(digits=6)       # significant digits

    # log.coa = results of PCA limited to the nbaxes first factorial axis
    datLog=log.coa$ind$coord
  }

  if (length(objects()[-which(objects() %in% c('methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
  } 



  if(pcaYesNo=="nopca"){
    datLog=datSpecies
    tabInertia="No PCA"
  }

  if (length(objects()[-which(objects() %in% c('methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
  }




# STEP 3 : CLUSTERING OF THE LOGEVENTS

print("######## STEP 3 CLUSTERING ########")

  if(methMetier=="hac"){
    # Distance between individuals
    log.dist <- dist(datLog, method=param3)             
    # HAC (Hierarchical Ascendant Classification)
    log.hac <- hclust(log.dist, method=param4)

    if (length(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
    }

    # Determine the number of cluster thanks to the scree-test
    inerties.vector=log.hac$height[order(log.hac$height,decreasing=T)]
    nb.finalclusters=which(scree(inerties.vector)$epsilon<0)[1]

    # Cut the dendogram at the selected level
    clusters=cutree(log.hac,k=nb.finalclusters)

    # Calculation of each cluster size
    n=dim(dat)[1]
    effgroupe=numeric()
    for(k in 1:nb.finalclusters){
      effgroupe[k]=length(which(clusters==k))
    }   
                                                
    # Rectangles plotting
    png("CAHdendogram.png", width = 1200, height = 800)
    plot(log.hac)
    rect.hclust(log.hac, k=nb.finalclusters)
    dev.off()

    resval=test.values(clusters,datSpecies)
    target=targetspecies(resval)


    # Mean profil of the dataset
    meanprofile=mean(datSpecies)
    png("Mean profile of the dataset.png", width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Mean profile of the dataset", xlab="Species", ylab="Euros")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    dev.off()

    # Mean profils by cluster
    nbClust=nb.finalclusters
    nbSpec=ncol(datSpecies)
    mprofil=numeric()
    for(i in 1:nbClust){
      mprofilclusti=mean(datSpecies[which(clusters==i),])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png("Mean profile by cluster.png", width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(mprofil[i,], main=paste("Cluster",i), xlab="Species", ylab="Euros")
      par(op2)
    }
    par(op)
    title(main="Mean profile by cluster")
    dev.off()

    # Number of Logevents by cluster
    x=c(1:nbClust)
    png("Number of Logevents by cluster.png", width = 1200, height = 800)
    barplot(effgroupe, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    dev.off()

    # Target Species profils (test-value)
    targetresval=numeric()
    nameTargetPlot=character()
    for(i in 1:nbClust){
      nomtargeti=as.character(target$tabnomespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      numtargeti=as.numeric(target$tabnumespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      nameTargetPloti=rep("",nbSpec)
      nameTargetPloti[numtargeti]=nomtargeti
      nameTargetPlot=rbind(nameTargetPlot,nameTargetPloti)
      targetresvalclusti=rep(0,nbSpec)
      targetresvalclusti[numtargeti]=resval[nomtargeti,i]
      targetresval=rbind(targetresval,targetresvalclusti)
    }

    png("Profile of target species by cluster.png", width = 1200, height = 800)
    #op <- par(mfrow=c(rep(ceiling(sqrt(nbClust)),2)))
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,],names.arg=nameTargetPlot[i,], main=paste("Cluster",i), xlab="Species", ylab="Test-value")
      par(op2)
    }
    par(op)
    title(main="Profile of target species by cluster")
    dev.off()


    return(list(clusters=clusters, effgroupe=effgroupe, datSpecies=datSpecies, tabInertia=tabInertia, datLog=datLog, nameTarget=target$tabnomespcib))

  }

  if (length(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
  }



  if(methMetier=="kmeans"){
    # Calculation of optimal k thanks to the variance intraclass
    varintra <- numeric()
    for (k in 2:10){
      clustersKmeans<-kmeans(datLog, k, iter.max=20, nstart=5)
      varintra [k] <- sum(clustersKmeans$withinss)
    }
    png("varintrakmeans.png", width = 1200, height = 800)
    plot(varintra)
    dev.off()

    if (length(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
    }

    diffvarintra=diff(varintra,na.rm=T)
    diffdiffvar=diff(diffvarintra,na.rm=T)
    maxdiffdiff=max(diffdiffvar,na.rm=T)
    k=which(diffdiffvar==maxdiffdiff)+1

    # KMEANS
    clusters<-kmeans(datLog, k, iter.max=20, nstart=6, algorithm="Hartigan-Wong")

    resval=test.values(clusters$cluster,datSpecies)
    target=targetspecies(resval)
    
    #### Inertie
    
    
    # Mean profil of the dataset
    meanprofile=mean(datSpecies)
    png("Mean profile of the dataset.png", width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Mean profile of the dataset", xlab="Species", ylab="Euros")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    dev.off()

    # Mean profils by cluster
    nbClust=length(clusters$size)
    nbSpec=ncol(datSpecies)
    mprofil=numeric()
    for(i in 1:nbClust){
      mprofilclusti=mean(datSpecies[which(clusters$cluster==i),])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png("Mean profile by cluster.png", width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(mprofil[i,], main=paste("Cluster",i), xlab="Species", ylab="Euros")
      par(op2)
    }
    par(op)
    title(main="Mean profile by cluster")
    dev.off()

    # Number of Logevents by cluster
    x=c(1:nbClust)
    png("Number of Logevents by cluster.png", width = 1200, height = 800)
    barplot(clusters$size, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    dev.off()

    # Target Species profils (test-value)
    targetresval=numeric()
    nameTargetPlot=character()
    for(i in 1:nbClust){
      nomtargeti=as.character(target$tabnomespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      numtargeti=as.numeric(target$tabnumespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      nameTargetPloti=rep("",nbSpec)
      nameTargetPloti[numtargeti]=nomtargeti
      nameTargetPlot=rbind(nameTargetPlot,nameTargetPloti)
      targetresvalclusti=rep(0,nbSpec)
      targetresvalclusti[numtargeti]=resval[nomtargeti,i]
      targetresval=rbind(targetresval,targetresvalclusti)
    }

    png("Profile of target species by cluster.png", width = 1200, height = 800)
    #op <- par(mfrow=c(rep(ceiling(sqrt(nbClust)),2)))
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,],names.arg=nameTargetPlot[i,], main=paste("Cluster",i), xlab="Species", ylab="Test-value")
      par(op2)
    }
    par(op)
    title(main="Profile of target species by cluster")
    dev.off()

    
    return(list(clusters=clusters, datSpecies=datSpecies, tabInertia=tabInertia, datLog=datLog, nameTarget=target$tabnomespcib))

  }                                                                                                         

  if (length(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
  }



  if(methMetier=="pam"){

    # Calculation of optimal k thanks to the silhouette
    clustersPam.silcoeff <- numeric()
    for (k in 2:10){
      clustersPam=pam(datLog,k)
      clustersPam.silcoeff [k] <- clustersPam$silinfo$avg.width
    }

    png("silcoeffpam.png", width = 1200, height = 800)
    plot(clustersPam.silcoeff)     # k optimal corresponds to maximum of silhouette's coefficients
    dev.off()

    if (length(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
    }

    clustersPam.silcoeff
    max=max(clustersPam.silcoeff, na.rm=T)
    k=which(clustersPam.silcoeff==max)
    clusters=pam(datLog,k)   # PAM with optimal k
    summary(clusters)

    # Silhouette of PAM with optimal k
    clusters.si=silhouette(clustersPam,do.col.sort)

    plot(clusters.si, col=rainbow(length(clusters$id.med)), main="Silhouette of Pam's clustering")
    savePlot(filename="Silhouette of Pam's Clustering", type='png', restoreConsole = TRUE)
    dev.off()

    resval=test.values(clusters$clustering,datSpecies)
    target=targetspecies(resval)
    

    # Mean profil of the dataset
    meanprofile=mean(datSpecies)
    png("Mean profile of the dataset.png", width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Mean profile of the dataset", xlab="Species", ylab="Euros")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    dev.off()

    # Mean profils by cluster
    nbClust=length(clusters$id.med)
    nbSpec=ncol(datSpecies)
    mprofil=numeric()
    for(i in 1:nbClust){
      mprofilclusti=mean(datSpecies[which(clusters$clustering==i),])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png("Mean profile by cluster.png", width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(mprofil[i,], main=paste("Cluster",i), xlab="Species", ylab="Euros")
      par(op2)
    }
    par(op)
    title(main="Mean profile by cluster")
    dev.off()

    # Number of Logevents by cluster
    x=c(1:nbClust)
    png("Number of Logevents by cluster.png", width = 1200, height = 800)
    barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    dev.off()

    # Target Species profils (test-value)
    targetresval=numeric()
    nameTargetPlot=character()
    for(i in 1:nbClust){
      nomtargeti=as.character(target$tabnomespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      numtargeti=as.numeric(target$tabnumespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      nameTargetPloti=rep("",nbSpec)
      nameTargetPloti[numtargeti]=nomtargeti
      nameTargetPlot=rbind(nameTargetPlot,nameTargetPloti)
      targetresvalclusti=rep(0,nbSpec)
      targetresvalclusti[numtargeti]=resval[nomtargeti,i]
      targetresval=rbind(targetresval,targetresvalclusti)
    }

    png("Profile of target species by cluster.png", width = 1200, height = 800)
    #op <- par(mfrow=c(rep(ceiling(sqrt(nbClust)),2)))
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,],names.arg=nameTargetPlot[i,], main=paste("Cluster",i), xlab="Species", ylab="Test-value")
      par(op2)
    }
    par(op)
    title(main="Profile of target species by cluster")
    dev.off()

    return(list(clusters=clusters, datSpecies=datSpecies, tabInertia=tabInertia, datLog=datLog, nameTarget=target$tabnomespcib))

  }

  if (length(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
  }



  if(methMetier=="clara"){
    n=dim(datLog)[1]

    # Calculation of optimal k thanks to the silhouette
    clustersClara.silcoeff <- numeric()
    for (k in 2:10){
      clustersClara=clara(datLog, k, metric=param3, stand=F, samples=5, sampsize=min(n,40+2*k))
      clustersClara.silcoeff [[k]] <- clustersClara$silinfo$avg.width
    }

    png("silcoeffclara.png", width = 1200, height = 800)
    plot(clustersClara.silcoeff)               # k optimal corresponds to maximum of silhouette's coefficients
    dev.off()
    
    if (length(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])   
    }

    clustersClara.silcoeff
    max=max(clustersClara.silcoeff, na.rm=T)
    k=which(clustersClara.silcoeff==max)
    clusters=clara(datLog, k, metric=param3, stand=F, samples=5, sampsize=min(n,40+2*k))  # CLARA with optimal k
    summary(clusters)

    #Silhouette of CLARA with optimal k
    clusters.si=silhouette(clusters, full=T)
    
    plot(clusters.si, col=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Silhouette of Clara's clustering obtained with the best sample")  # full silhouette obtained from the best sample
    savePlot(filename="Silhouette of Clara's Clustering", type='png', restoreConsole = TRUE)
    dev.off()

# Graph cluster
#  clusplot(clustersClara,col.p = clustersClara$clustering)

    resval=test.values(clusters$cluster,datSpecies)
    target=targetspecies(resval)
    
    
    # Mean profil of the dataset
    meanprofile=mean(datSpecies)
    png("Mean profile of the dataset.png", width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Mean profile of the dataset", xlab="Species", ylab="Euros")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    dev.off()

    # Mean profils by cluster
    nbClust=length(clusters$i.med)
    nbSpec=ncol(datSpecies)
    mprofil=numeric()
    for(i in 1:nbClust){
      mprofilclusti=mean(datSpecies[which(clusters$clustering==i),])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png("Mean profile by cluster.png", width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(mprofil[i,], main=paste("Cluster",i), xlab="Species", ylab="Euros")
      par(op2)
    }
    par(op)
    title(main="Mean profile by cluster")
    dev.off()

    # Number of Logevents by cluster
    x=c(1:nbClust)
    png("Number of Logevents by cluster.png", width = 1200, height = 800)
    barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    dev.off()

    # Target Species profils (test-value)
    targetresval=numeric()
    nameTargetPlot=character()
    for(i in 1:nbClust){
      nomtargeti=as.character(target$tabnomespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      numtargeti=as.numeric(target$tabnumespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      nameTargetPloti=rep("",nbSpec)
      nameTargetPloti[numtargeti]=nomtargeti
      nameTargetPlot=rbind(nameTargetPlot,nameTargetPloti)
      targetresvalclusti=rep(0,nbSpec)
      targetresvalclusti[numtargeti]=resval[nomtargeti,i]
      targetresval=rbind(targetresval,targetresvalclusti)
    }

    png("Profile of target species by cluster.png", width = 1200, height = 800)
    #op <- par(mfrow=c(rep(ceiling(sqrt(nbClust)),2)))
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,],names.arg=nameTargetPlot[i,], main=paste("Cluster",i), xlab="Species", ylab="Test-value")
      par(op2)
    }
    par(op)
    title(main="Profile of target species by cluster")
    dev.off()
    
    return(list(clusters=clusters, datSpecies=datSpecies, tabInertia=tabInertia, datLog=datLog, nameTarget=target$tabnomespcib))

  }
  
  if (length(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])  
  }
  # end of the methods

} # end of the function "classif"

