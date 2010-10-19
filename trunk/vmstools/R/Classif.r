##############################################################################################################
#                                                                                                            #
#     CLASSIF : FUNCTION FOR SPECIES SELECTION, PCA, CLASSIFICATION OF LOGEVENTS AND MÉTIERS COMPUTATION     #
#                                                                                                            #
##############################################################################################################

                                                                                
#############################################################
#  PARAMETERS TO FILL IN DEPENDING OF THE SELECTED METHODS  #
#############################################################
#
# STEP 1 : Choose the method for the selection of species : 3 methods : HAC, totale, logevent
#
# HAC : we view as main species, all species which are not belonging to the cluster with the smaller mean
# 2 parameters to fill in : param1 for the distance between individuals ("euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski")
#                           param2 for the agregative distance (distance between clusters) ("ward", "single", "complete", "average", "mcquitty", "median" or "centroid")
#
# totale : we take all species by decreasing importance order until having a given percentage (70, 80, 90%) of the total capture
# 1 parameter to fill in : param1 = the threshold (70, 80, 90...)
#
# logevent : we take all species which represents at least a given percentage (10, 20, 30%) for at least a logevent
# 1 parameter to fill in : param1 = the threshold (10, 20, 30...)
#
#
# STEP 2 : PCA or not ?
#
# pca : we execute a PCA on the individuals "logevents", and compute the classification on the factorial coordinates of the logevents
# 1 parameter to fill in : criterion for the criterion of number of factorial axis selection (70percents for taking all first factorial axis until having at least 70% of initial inertia)
#                                                                                            (screetest for using the scree-test)
#
# nopca : we choose to work with the raw data to compute the classification of logevents
#
#
# STEP 3 : Choose the method for the classification of logevents : HAC, KMEANS, PAM, CLARA ?
#
# HAC, 2 parameters to fill in : param3 for the distance between individuals ("euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski")
#                                param4 for the agregative distance (distance between clusters) ("ward", "single", "complete", "average", "mcquitty", "median" ou "centroid")
#
# KMEANS, no parameter to fill in
#
# PAM, no parameter to fill in
#
# CLARA, 1 parameter to fill in : param3 for the dissimilarity between individuals ("euclidean" or "manhattan")




# Classif : function for the selection of species (hac, totale, logevent), PCA (pca, nopca), classification of logevents (hac, kmeans, pam, and clara) and métiers computation (thanks to test-values)
classif=function(dat,analysisName="",methSpecies="hac",param1="euclidean",param2="ward",pcaYesNo="pca",criterion="70percents",methMetier="clara",param3="euclidean",param4=NULL){


# First simplify the names of columns
names(dat)[-1] <- unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))
names(dat)




##############################################################################################################################################
#                                                STEP 1 : SELECTION OF SPECIES                                                               #
##############################################################################################################################################


print("######## STEP 1 SELECTION OF MAIN SPECIES ########")

  if(methSpecies=="hac"){

    toutfait=FALSE
    p=ncol(dat)   # Number of species
 
    # Transform quantities to proportions of total quantity caught by logevent
    propdat=transformation_proportion(dat[,2:p])
    names(propdat)
    
    # Transposing data
    table_var=table_variables(propdat)

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
#    CAH
#    esp_dist=dist(table_var, method=param1)
#    cah_var=hclust(esp_dist, method=param2)

    # HAC
    cah_var=hcluster(table_var, method=param1, link=param2)

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
    # Select the number of clusters by scree-test
    inerties.vector=cah_var$height[order(cah_var$height,decreasing=T)]
    nb.finalclusters=which(scree(inerties.vector)$epsilon<0)[1]

    # Dendogram cutting at the selected level
    cah_cluster_var=cutree(cah_var,k=nb.finalclusters)

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
    # Selection of main species
    espprinc=select_species(dat[,2:p],cah_cluster_var)
    espdata=names(dat[,2:p])
    
    # Selection of possible main species among residuals species
    while(toutfait==FALSE){
      ind_princ=which(is.element(espdata,espprinc))
      ind_autres=setdiff(espdata,espprinc)
      princ=table_var[ind_princ,]
      resi=table_var[ind_autres,]

      # CAH on residual species group
#      resiesp_dist=dist(resi, method=param1)
#      resicah_var=hclust(resiesp_dist, method=param2)

      # HAC on residual species group
      resicah_var=hcluster(resi, method=param1, link=param2)
      
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
#      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
#      gc()
    }

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
    # Return only the main species of the dataset
    datSpecies=building_tab_pca(propdat,espprinc)
    #write.table(datSpecies, file="datSpecies.txt", quote=T, dec='.', sep=';', col.names=T, row.names=F)
    
  }
  
  
  if(methSpecies=="totale"){
  
    p=length(dat[1,])   # Number of species
    propdat=transformation_proportion(dat[,2:p])
    
    # Total quantity caught species by species
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
    
    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
    # We are taking all species until having at least param1% of total catch
    seuil=param1
    pourcent=which(propespcum<seuil)
    espsel=numesp[1:(length(pourcent)+1)]
    # We are taking the name of selected species
    nomespsel=names(dat[espsel])
    # We are bulding the table with main species and aggregated other species
    datSpecies=building_tab_pca(propdat,nomespsel)

  }


  if(methSpecies=="logevent"){

    p=length(dat[1,])   # Number of species
    # Transform quantities to proportions of total quantity caught by logevent
    propdat=transformation_proportion(dat[,2:p])
    seuil=param1
    # Selection of species making up over param1% of logevent's captures
    pourcent <- apply(propdat,1,function(x) which(x>seuil))   
    nomespsel <- names(propdat)[unique(unlist(pourcent))]  
    
    # We are bulding the table with main species and aggregated other species
    datSpecies=building_tab_pca(propdat,nomespsel)

  }







##############################################################################################################################################
#                                                          STEP 2 : PCA OR NOT                                                               #
##############################################################################################################################################



print("######## STEP 2 PCA/NO PCA ON CATCH PROFILES ########")

  if(pcaYesNo=="pca"){
    # PCA (Principal Component Analysis)
    log.pca <- PCA(datSpecies, graph=T)
    
    savePlot(filename=paste(analysisName,'Species projection on the two first factorial axis',sep="_"), type='png', restoreConsole = TRUE)
    dev.off()
    savePlot(filename=paste(analysisName,'Individuals projection on the two first factorial axis',sep="_"), type='png', restoreConsole = TRUE)
    dev.off()

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
    # Data frame given eigenvalues, inertia and cumulative inertia of factorial axis
    tabInertia=data.frame(cbind(Axis=1:length(log.pca$eig[,1]), Eigenvalues=log.pca$eig[,1], Inertia=log.pca$eig[,2], CumulativeInertia=log.pca$eig[,3]))                  

    # Determine the number of axis to keep
    if(criterion=="70percents"){
      nbaxes=which(log.pca$eig[,3]>70)[1]   # we are taking the axis until having 70% of total inertia
    }
    # OR
    if(criterion=="screetest"){
      nbaxes=which(scree(tabInertia[,3])$epsilon<0)[1]  # thanks to the scree-test
    }
      
    # Eigenvalues and relative graphics
    log.pca$eig                                                  

    png(paste(analysisName,"Eigen values.png",sep="_"), width = 1200, height = 800)
    x=1:length(log.pca$eig[,1])
    barplot(log.pca$eig[,1],names.arg=x, main="Eigen values")
    dev.off()
    png(paste(analysisName,"Percentage of Inertia.png",sep="_"), width = 1200, height = 800)
    color=rep("grey",length(log.pca$eig[,1]))
    color[1:nbaxes]="green"
    barplot(log.pca$eig[,2],names.arg=x, col=color, main="Percentage of Inertia of factorial axis", xlab="Axis", ylab="% of Inertia")
    dev.off()
    png(paste(analysisName,"Cumulative Percentage of Inertia.png",sep="_"), width = 1200, height = 800)
    color=rep("grey",length(log.pca$eig[,1]))
    numaxe=min(which(log.pca$eig[,3]>70))
    color[1:numaxe]="green"
    barplot(log.pca$eig[,3],names.arg=x, col=color, main="Cumulative Percentage of Inertia of factorial axis", xlab="Axis", ylab="% of Inertia")
    abline(h=70, col="red")
    text(1,72, "70% of Inertia", col = "red", adj = c(0, -.1))
    dev.off()   
    
    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
    # PCA with the good number of axis
    log.coa=PCA(datSpecies, graph=F, ncp=nbaxes)
    options(digits=6)       # significant digits
    
    # Projection of variables Species on the first factorial axis
    png(paste(analysisName,"Projection of Species on first factorial axis.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(log.pca,choix="var",axes = c(1, 2),new.plot=FALSE,lim.cos2.var = 0.3)
    plot(log.pca,choix="var",axes = c(2, 3),new.plot=FALSE,lim.cos2.var = 0.3)
    plot(log.pca,choix="var",axes = c(1, 3),new.plot=FALSE,lim.cos2.var = 0.3)
    plot(log.pca,choix="var",axes = c(1, 4),new.plot=FALSE,lim.cos2.var = 0.3)
    plot(log.pca,choix="var",axes = c(2, 4),new.plot=FALSE,lim.cos2.var = 0.3)
    plot(log.pca,choix="var",axes = c(3, 4),new.plot=FALSE,lim.cos2.var = 0.3)
    par(op)
    dev.off()

    # log.coa = results of PCA limited to the nbaxes first factorial axis
    datLog=signif(log.coa$ind$coord, 5)
    #write.table(datLog, file="datLog.txt", quote=T, dec='.', sep=';', col.names=T, row.names=F)
  }



  if(pcaYesNo=="nopca"){
    datLog=datSpecies
    tabInertia="No PCA"
  }







##############################################################################################################################################
#                                               STEP 3 : CLUSTERING OF THE LOGEVENTS                                                         #
##############################################################################################################################################



########################################################################################################################################   HAC

print("######## STEP 3 CLUSTERING ########")

  if(methMetier=="hac"){
  
    classifWithinVar=numeric()
    classifBetweenVar=numeric()
    classifQuality=numeric()
    sampleList=numeric()
    mProfilSample=numeric()

    #totaleVarPart=numeric()   #(used for the calcul of totale variance)
    nbLog=nrow(datLog)
    nbDim=ncol(datLog)
    
    # Center of gravity of datLog
    centerOfGravityDatLog=colMeans(datLog)
    
    # HAC like CLARA (HAC on sample, affectation of each logevent to a cluster, quality of classification, do it 5 times, choose the sample which gives the best quality of classification)
    for(i in 1:5){
      # Sample of size 15000 logevents
      sam=sample(1:nbLog,size=15000,replace=F)
      # Record the 5 samples
      sampleList=rbind(sampleList,sam)
      outofsam=setdiff(1:nbLog,sam)
      sampleDatLog=datLog[sam,]

      # HAC on the sample
      log.hac=hcluster(sampleDatLog, method=param3, link=param4)
      inerties.vector=log.hac$height[order(log.hac$height,decreasing=T)]
      nbClust=which(scree(inerties.vector)$epsilon<0)[1]
    
      # Cut the dendogram at the selected level
      sampleClusters=cutree(log.hac,k=nbClust)

      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
      gc()

      # Add the cluster to each logevent of the sample
      sampleDatLogWithClusters=cbind(sampleDatLog,sampleClusters)

      # Discriminante analysis on the learning dataset "sampleDatLogWithClusters"
      learning.lda=lda(sampleDatLogWithClusters[,1:nbDim],sampleDatLogWithClusters[,ncol(sampleDatLogWithClusters)])

      # Other Logevents (which aren't belonging to the sample)
      otherLog=datLog[outofsam,]

      # Predict the cluster for each logevent of "otherLog" thanks to a discriminante analysis
      result.lda=predict(learning.lda,otherLog)

      # Add the cluster to each logevent of "otherLog"
      otherDatLogWithClusters=cbind(otherLog, result.lda$class)
      colnames(otherDatLogWithClusters)=colnames(sampleDatLogWithClusters)

      # Rebuilt datLog with clusters
      datLogWithClusters=rbind(sampleDatLogWithClusters,otherDatLogWithClusters)

      # Clusters for logevents
      clusters=datLogWithClusters[,ncol(datLogWithClusters)]

      # Between variance of classification
      sizeClusters=numeric()
      centerOfGravityClassif=numeric()
      centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)
      for(k in 1:nbClust){
        clusti=datLogWithClusters[which(clusters==k),1:nbDim]
        sizeClusters[k]=nrow(clusti)
        centerOfGravityClusti=colMeans(clusti)
        centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
      }
      classifBetweenVar=cbind(classifBetweenVar,1/nbLog*sum(sizeClusters*((dist(centerOfGravityClassif)[1:nbClust])^2)))

      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
      gc()
    
      # Within variance of classification
      withinVarClusters=numeric()
      for(k in 1:nbClust){  # Within variance by cluster

        clusti=datLogWithClusters[which(clusters==k),1:nbDim]
        sizeClusti=nrow(clusti)
        centerOfGravityClusti=colMeans(clusti)
        withinVarClustiPart=numeric()
        withinVarClusters[k]=1/sizeClusti*sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))

      }
      classifWithinVar=cbind(classifWithinVar,1/nbLog*sum(sizeClusters*withinVarClusters))   # Within variance of classification
      
      # Mean profiles by cluster for each sample
      nbSpec=ncol(datSpecies)
      mprofil=numeric()
      blank=rep(00000000000,nbSpec)
      for(k in 1:nbClust){
        mprofilclusti=mean(datSpecies[which(clusters==k),])
        mprofil=rbind(mprofil,mprofilclusti)
      }
      mprofil=rbind(mprofil,blank)
      
      mProfilSample=rbind(mProfilSample,mprofil)

      Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
      gc()
      
    }

    
    # Classification's quality
    classifQuality=classifWithinVar/classifBetweenVar
    which.min(classifQuality)

    # Select the sample which gives the smaller classification's quality (the best sample)
    sam=sampleList[which.min(classifQuality),]
    outofsam=setdiff(1:nbLog,sam)
    sampleDatLog=datLog[sam,]

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
    # HAC with the best sample
    log.hac=hcluster(sampleDatLog, method=param3, link=param4)
    
    # Determine the number of cluster thanks to the scree-test
    inerties.vector=log.hac$height[order(log.hac$height,decreasing=T)]
    nbClust=which(scree(inerties.vector)$epsilon<0)[1]

    # Cut the dendogram at the selected level
    sampleClusters=cutree(log.hac,k=nbClust)
    
    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()

    sampleDatLogWithClusters=cbind(sampleDatLog,sampleClusters)

    nbLogSample=nrow(sampleDatLog)
    nbDim=ncol(sampleDatLog)

    # Discriminante analysis on the learning dataset
    learning.lda=lda(sampleDatLogWithClusters[,1:nbDim],sampleDatLogWithClusters[,ncol(sampleDatLogWithClusters)])

    otherLog=datLog[outofsam,]

    # Predict the cluster for the other logevent
    result.lda=predict(learning.lda,otherLog)

    otherDatLogWithClusters=cbind(otherLog, result.lda$class)
    colnames(otherDatLogWithClusters)=colnames(sampleDatLogWithClusters)

    # Rebuilt complete datLog with clusters
    datLogWithClusters=rbind(sampleDatLogWithClusters,otherDatLogWithClusters)
    clusters=datLogWithClusters[,ncol(datLogWithClusters)]

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
    # Within and between variance of clusters and classification
    centerOfGravityClassif=numeric()
    withinVarClusters=numeric()
    centerOfGravityDatLog=colMeans(datLog)
    centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)
    for(k in 1:nbClust){  # Within variance by cluster

      clusti=datLogWithClusters[which(clusters==k),1:nbDim]
      centerOfGravityClusti=colMeans(clusti)
      centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
      sizeClusti=nrow(clusti)
      centerOfGravityClusti=colMeans(clusti)
      withinVarClustiPart=numeric()
      withinVarClusters[k]=1/sizeClusti*sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))

    }
    # Between variance
    classifBetweenVar=1/nbLog*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
    # Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
    withinVarClusterOnTot=withinVarClusters/(classifBetweenVar+sum(withinVarClusters))*100
    betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+sum(withinVarClusters))*100


    # Calculation of each cluster size
    n=nrow(datLog)
    sizeClusters=numeric()
    for(k in 1:nbClust){
      sizeClusters[k]=length(which(clusters==k))
    }

    # Compute the test-values for species
    resval=test.values(clusters,datSpecies)
    # Determine the target species
    target=targetspecies(resval)

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()


    # Projections on the first factorial plans
    png(paste(analysisName,"HAC_Projections.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
    plot(datLog[,1], datLog[,3], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
    plot(datLog[,1], datLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
    plot(datLog[,2], datLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
    plot(datLog[,3], datLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
    par(op)
    dev.off()

    
    # Rectangles plotting
    png(paste(analysisName,"HAC_Dendogram.png",sep="_"), width = 1200, height = 800)
    plclust(log.hac,labels=F,hang=-1,ann=F)
    title(main="HAC dendogram",xlab="Logevents",ylab="Height")
    rect.hclust(log.hac, k=nbClust)
    dev.off()


    # Mean profile of the dataset
    meanprofile=mean(datSpecies)
    png(paste(analysisName,"Mean profile of the dataset.png",sep="_"), width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Mean profile of the dataset", xlab="Species", ylab="Percentage of catch")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    dev.off()
    

    # Mean profiles by cluster
    nbSpec=ncol(datSpecies)
    mprofil=numeric()
    for(k in 1:nbClust){
      mprofilclusti=mean(datSpecies[which(clusters==k),])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png(paste(analysisName,"Mean profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(mprofil[i,], cex.names=1, xlab="Species", ylab="Percentage of catch")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main="Mean profile by cluster")
    dev.off()
    
    
    # Standard deviation profile by cluster
    sdprofil=numeric()
    for(i in 1:nbClust){
      sdprofilclusti=sd(datSpecies[which(clusters==i),])
      sdprofil=rbind(sdprofil,sdprofilclusti)
    }
    png(paste(analysisName,"Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(sdprofil[i,], cex.names=1, xlab="Species", ylab="Percentage of catch")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main="Standard deviation profile by cluster")
    dev.off()

#    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
#    gc()
    
    
    # Number of Logevents by cluster
    x=c(1:nbClust)
    png(paste(analysisName,"Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
    coord=barplot(sizeClusters, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    barplot(sizeClusters, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
    text(coord,sizeClusters-500,sizeClusters,font=2)
    text(coord,sizeClusters+500,sizeClusters,font=2)
    dev.off()


    # Target Species profiles (test-value)
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

    png(paste(analysisName,"Profile of target species by cluster.png",sep="_"), width = 1200, height = 800)
    #op <- par(mfrow=c(rep(ceiling(sqrt(nbClust)),2)))
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,],names.arg=nameTargetPlot[i,], cex.names=1, xlab="Species", ylab="Test-value")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main="Profile of target species by cluster")
    dev.off()
    
    return(list(clusters=clusters, sizeClusters=sizeClusters, datSpecies=datSpecies, tabInertia=tabInertia, datLog=datLog, nameTarget=target$tabnomespcib, betweenVarClassifOnTot=betweenVarClassifOnTot, mProfilSample=mProfilSample))

  }





########################################################################################################################################   KMEANS

  if(methMetier=="kmeans"){
    # Calculation of optimal k thanks to within variance
    varintra=numeric()
    nbLog=nrow(datLog)
    for (k in 2:15){
      clustersKmeans=kmeans(datLog, k, iter.max=20, nstart=5)
      varintra[k]=1/nbLog*sum(clustersKmeans$withinss)
    }
    png(paste(analysisName,"Within Variance Kmeans.png",sep="_"), width = 1200, height = 800)
    plot(varintra,main="Within clusters variance",xlab="Number of clusters",ylab="Within Variance")
    dev.off()

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()

    diffvarintra=diff(varintra,na.rm=T)
    diffdiffvar=diff(diffvarintra,na.rm=T)
    maxdiffdiff=max(diffdiffvar,na.rm=T)
    k=which(diffdiffvar==maxdiffdiff)+1

    # KMEANS with k optimal
    clusters<-kmeans(datLog, k, iter.max=20, nstart=6, algorithm="Hartigan-Wong")


    # Quality of classification
    centerOfGravityClassif=numeric()
    centerOfGravityDatLog=colMeans(datLog)
    centerOfGravityClassif=rbind(centerOfGravityDatLog,clusters$centers)
    classifBetweenVar=1/nbLog*sum(clusters$size*((dist(centerOfGravityClassif)[1:k])^2))
    classifWithinVar=1/nbLog*sum(clusters$withinss)
    classifQuality=classifWithinVar/classifBetweenVar
    withinVarClustersOnTot=(clusters$withinss/nbLog)*100/(classifBetweenVar+classifWithinVar)
    betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+classifWithinVar)*100

    # Compute the test-values for species
    resval=test.values(clusters$cluster,datSpecies)
    #Determine the target species
    target=targetspecies(resval)

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()


    # Projections on the first factorial plans
    png(paste(analysisName,"KMEANS_Projections.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
    plot(datLog[,1], datLog[,3], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
    plot(datLog[,1], datLog[,4], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
    plot(datLog[,2], datLog[,4], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
    plot(datLog[,3], datLog[,4], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
    par(op)
    dev.off()
    
    
    # Mean profile of the dataset
    meanprofile=mean(datSpecies)
    png(paste(analysisName,"Mean profile of the dataset.png",sep="_"), width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Mean profile of the dataset", xlab="Species", ylab="Percentage of catch")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    dev.off()


    # Mean profiles by cluster
    nbClust=length(clusters$size)
    nbSpec=ncol(datSpecies)
    mprofil=numeric()
    for(i in 1:nbClust){
      mprofilclusti=mean(datSpecies[which(clusters$cluster==i),])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png(paste(analysisName,"Mean profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(mprofil[i,], cex.names=1, main=paste("Cluster",i), xlab="Species", ylab="Percentage of catch")
      par(op2)
      #mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main="Mean profile by cluster")
    dev.off()
    
#    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
#    gc()
    
    
    # Standard deviation profiles by cluster
    sdprofil=numeric()
    for(i in 1:nbClust){
      sdprofilclusti=sd(datSpecies[which(clusters$cluster==i),])
      sdprofil=rbind(sdprofil,sdprofilclusti)
    }
    png(paste(analysisName,"Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(sdprofil[i,], cex.names=1, main=paste("Cluster",i), xlab="Species", ylab="Percentage of catch")
      par(op2)
      #mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main="Standard deviation profile by cluster")
    dev.off()


    # Number of Logevents by cluster
    x=c(1:nbClust)
    png(paste(analysisName,"Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
    coord=barplot(clusters$size, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    barplot(clusters$size, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
    text(coord,clusters$size-600,clusters$size,font=2)
    text(coord,clusters$size+600,clusters$size,font=2)
    dev.off()



    # Target Species profiles (test-value)
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

    png(paste(analysisName,"Profile of target species by cluster.png",sep="_"), width = 1200, height = 800)
    #op <- par(mfrow=c(rep(ceiling(sqrt(nbClust)),2)))
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,], cex.names=1, main=paste("Cluster",i), names.arg=nameTargetPlot[i,], xlab="Species", ylab="Test-value")
      par(op2)
      #mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main="Profile of target species by cluster")
    dev.off()

    
    return(list(clusters=clusters, datSpecies=datSpecies, tabInertia=tabInertia, datLog=datLog, nameTarget=target$tabnomespcib, betweenVarClassifOnTot=betweenVarClassifOnTot))

  }                                                                                                         





########################################################################################################################################   PAM

  if(methMetier=="pam"){

    # Calculation of optimal k thanks to the silhouette
    clustersPam.silcoeff <- numeric()
    for (k in 2:10){
      clustersPam=pam(datLog,k)
      clustersPam.silcoeff [k] <- clustersPam$silinfo$avg.width
    }

    png(paste(analysisName,"silcoeffpam.png",sep="_"), width = 1200, height = 800)
    plot(clustersPam.silcoeff)     # k optimal corresponds to maximum of silhouette's coefficients
    dev.off()

    clustersPam.silcoeff
    max=max(clustersPam.silcoeff, na.rm=T)
    k=which(clustersPam.silcoeff==max)
    
    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
    # PAM with optimal k
    clusters=pam(datLog,k)
    summary(clusters)

    # Silhouette of PAM with optimal k
#    clusters.si=silhouette(clustersPam,do.col.sort)
#
#    plot(clusters.si, col=rainbow(length(clusters$id.med)), main="Silhouette of Pam's clustering")
#    savePlot(filename="Silhouette of Pam's Clustering", type='png', restoreConsole = TRUE)
#    dev.off()

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()
    
    # Quality of classification
    nbLog=nrow(datLog)
    centerOfGravityClassif=numeric()
    centerOfGravityDatLog=colMeans(datLog)
    centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)

    for(i in 1:k){
      clusti=datLog[which(clusters$clustering==i),]
      centerOfGravityClusti=colMeans(clusti)
      centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
    }
    classifBetweenVar=1/nbLog*sum(clusters$clusinfo[,1]*((dist(centerOfGravityClassif)[1:k])^2))
    classifWithinVar=1/nbLog*sum(clusters$clusinfo[,1]*clusters$clusinfo[,3])
    classifQuality=classifWithinVar/classifBetweenVar
    withinVarClustersOnTot=(clusters$clusinfo[,1]*clusters$clusinfo[,3]/nbLog)*100/(classifWithinVar+classifBetweenVar)
    betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+classifWithinVar)*100
    

    # Compute the test-values for species
    resval=test.values(clusters$clustering,datSpecies)
    # Determine the target species
    target=targetspecies(resval)

    
    # Projections on the first factorial plans
    png(paste(analysisName,"PAM_Projections.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
    plot(datLog[,1], datLog[,3], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
    plot(datLog[,1], datLog[,4], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
    plot(datLog[,2], datLog[,4], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
    plot(datLog[,3], datLog[,4], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
    par(op)
    dev.off()

    
    # Mean profile of the dataset
    meanprofile=mean(datSpecies)
    png(paste(analysisName,"Mean profile of the dataset.png",sep="_"), width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Mean profile of the dataset", xlab="Species", ylab="Percentage of catch")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    dev.off()


    # Mean profile by cluster
    nbClust=length(clusters$id.med)
    nbSpec=ncol(datSpecies)
    mprofil=numeric()
    for(i in 1:nbClust){
      mprofilclusti=mean(datSpecies[which(clusters$clustering==i),])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png(paste(analysisName,"Mean profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(mprofil[i,], main=paste("Cluster",i), xlab="Species", ylab="Percentage of catch")
      par(op2)
    }
    par(op)
    title(main="Mean profile by cluster")
    dev.off()
    
#    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
#    gc()
    
    
    # Standard deviation profile by cluster
    sdprofil=numeric()
    for(i in 1:nbClust){
      sdprofilclusti=sd(datSpecies[which(clusters$clustering==i),])
      sdprofil=rbind(sdprofil,sdprofilclusti)
    }
    png(paste(analysisName,"Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(sdprofil[i,], xlab="Species", ylab="Percentage of catch")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main="Standard deviation profile by cluster")
    dev.off()


    # Number of Logevents by cluster
    x=c(1:nbClust)
    png(paste(analysisName,"Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
    coord=barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
    text(coord,clusters$clusinfo[,1]-600,clusters$clusinfo[,1],font=2)
    dev.off()


    # Target Species profile (test-value) by cluster
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

    png(paste(analysisName,"Profile of target species by cluster.png",sep="_"), width = 1200, height = 800)
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

    return(list(clusters=clusters, datSpecies=datSpecies, tabInertia=tabInertia, datLog=datLog, nameTarget=target$tabnomespcib, betweenVarClassifOnTot=betweenVarClassifOnTot))

  }





########################################################################################################################################   CLARA

  if(methMetier=="clara"){
    nbLog=nrow(datLog)

    # Calculation of optimal k thanks to the silhouette
    clustersClara.silcoeff=numeric()
    for (k in 2:15){
      clustersClara=clara(datLog, k, metric=param3, stand=F, samples=5, sampsize=min(nbLog,40+2*k))
      clustersClara.silcoeff[k]=clustersClara$silinfo$avg.width
    }

    png(paste(analysisName,"silcoeffclara.png",sep="_"), width = 1200, height = 800)
    plot(clustersClara.silcoeff, main="Silhouette of the classification", xlab="Number of clusters", ylab="Silhouette")               # k optimal corresponds to maximum of silhouette's coefficients
    dev.off()
    
    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()

    clustersClara.silcoeff
    max=max(clustersClara.silcoeff, na.rm=T)
    k=which(clustersClara.silcoeff==max)
    
    # CLARA with optimal k
    clusters=clara(datLog, k, metric=param3, stand=F, samples=5, sampsize=min(nbLog,40+2*k))  # CLARA with optimal k
    summary(clusters)

#    #Silhouette of CLARA with optimal k
#    clusters.si=silhouette(clusters, full=T)
#
#    plot(clusters.si, col=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Silhouette of Clara's clustering obtained with the best sample")  # full silhouette obtained from the best sample
#    savePlot(filename="Silhouette of Clara's Clustering", type='png', restoreConsole = TRUE)
#    dev.off()

# Graph cluster
#  clusplot(clustersClara,col.p = clustersClara$clustering)

    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc()


    # Quality of classification
    centerOfGravityClassif=numeric()
    centerOfGravityDatLog=colMeans(datLog)
    centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)

    for(i in 1:k){
      clusti=datLog[which(clusters$clustering==i),]
      centerOfGravityClusti=colMeans(clusti)
      centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
    }
    classifBetweenVar=1/nbLog*sum(clusters$clusinfo[,1]*((dist(centerOfGravityClassif)[1:k])^2))
    classifWithinVar=1/nbLog*sum(clusters$clusinfo[,1]*clusters$clusinfo[,3])
    classifQuality=classifWithinVar/classifBetweenVar
    withinVarClustersOnTot=(clusters$clusinfo[,1]*clusters$clusinfo[,3]/nbLog)*100/(classifWithinVar+classifBetweenVar)
    betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+classifWithinVar)*100


#    Store(objects()[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
#    gc()


    # Compute the test-values for species
    resval=test.values(clusters$cluster,datSpecies)
    # Determine the target species
    target=targetspecies(resval)
    
    # Projections on the first factorial plans
    png(paste(analysisName,"CLARA_Projections.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
    plot(datLog[,1], datLog[,3], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
    plot(datLog[,1], datLog[,4], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
    plot(datLog[,2], datLog[,4], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
    plot(datLog[,3], datLog[,4], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
    par(op)
    dev.off()


    # Mean profile of the dataset
    meanprofile=mean(datSpecies)
    png(paste(analysisName,"Mean profile of the dataset.png",sep="_"), width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Mean profile of the dataset", xlab="Species", ylab="Percentage of catch")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    dev.off()


    # Mean profile by cluster
    nbClust=length(clusters$i.med)
    nbSpec=ncol(datSpecies)
    mprofil=numeric()
    for(i in 1:nbClust){
      mprofilclusti=mean(datSpecies[which(clusters$clustering==i),])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png(paste(analysisName,"Mean profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(mprofil[i,], xlab="Species", ylab="Percentage of catch")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main="Mean profile by cluster")
    dev.off()
    
    
    # Standard deviation profile by cluster
    sdprofil=numeric()
    for(i in 1:nbClust){
      sdprofilclusti=sd(datSpecies[which(clusters$clustering==i),])
      sdprofil=rbind(sdprofil,sdprofilclusti)
    }
    png(paste(analysisName,"Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(sdprofil[i,], xlab="Species", ylab="Percentage of catch")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main="Standard deviation profile by cluster")
    dev.off()
    

    # Number of Logevents by cluster
    x=c(1:nbClust)
    png(paste(analysisName,"Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
    coord=barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
    text(coord,clusters$clusinfo[,1]-600,clusters$clusinfo[,1],font=2)
    text(coord,clusters$clusinfo[,1]+800,clusters$clusinfo[,1],font=2)
    dev.off()


    # Target Species profile (test-value) by cluster
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

    png(paste(analysisName,"Profile of target species by cluster.png",sep="_"), width = 1200, height = 800)
    #op <- par(mfrow=c(rep(ceiling(sqrt(nbClust)),2)))
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,],names.arg=nameTargetPlot[i,], xlab="Species", ylab="Test-value")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=F, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main="Profile of target species by cluster")
    dev.off()
    
    return(list(clusters=clusters, datSpecies=datSpecies, tabInertia=tabInertia, datLog=datLog, nameTarget=target$tabnomespcib, betweenVarClassifOnTot=betweenVarClassifOnTot))

  }  
  # end of the methods

} # end of the function "classif"