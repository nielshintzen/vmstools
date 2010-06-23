##############################################################################################################
#                                                                                                            #
#     CLASSIF : FUNCTION FOR SPECIES SELECTION, PCA, CLASSIFICATION OF LOGEVENTS AND MÉTIERS COMPUTATION     #
#                                                                                                            #
##############################################################################################################

                                                                                
#############################################################
#  PARAMETERS TO FILL IN DEPENDING OF THE SELECTED METHODS  #
#############################################################
#
# STEP 1 : Choose the method for the selection of species : 3 methods : CAH, totale, logevent
#
# CAH : we view as main species, all species which are not belonging to the cluster with the smaller mean
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
# acp : we execute a PCA on the individuals "logevents", and compute the classification on the factorial coordinates of the logevents
# noacp : we choose to work with the raw data to compute the classification of logevents
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




# Classif : function for the selection of species (cah, totale, logevent), PCA (acp, noacp), classification of logevents (cah, kmeans, pam, and clara) and métiers computation (thanks to test-values)
classif=function(dat,methspecies="cah",param1="euclidean",param2="ward",acpyesno="acp",methmetier="clara",param3="euclidean",param4=NULL){

#first siplify the names of columns 
names(dat)[-1] <- unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))
names(dat)

# STEP 1 : SELECTION OF SPECIES

print("######## STEP 1 SELECTION OF MAIN SPECIES ########")

  if(methspecies=="cah"){

    toutfait=FALSE
    p=length(dat[1,])   # Number of species
    # Transform quantities to proportions of total quantity caught by logevent
    InputData_profil=transformation_proportion(dat[,2:p])
    names(InputData_profil)
    # Transposing data
    table_var=table_variables(InputData_profil)

    Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   

    # CAH
    esp_dist=dist(table_var, method=param1)
    cah_var=hclust(esp_dist, method=param2)

    Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   

    # Select the number of clusters by scree-test
    inerties.vector=cah_var$height[order(cah_var$height,decreasing=T)]
    nb.finalclusters=which(scree(inerties.vector)$epsilon<0)[1]

    # Dendogram cutting at the selected level
    cah_cluster_var=cutree(cah_var,k=nb.finalclusters)

    Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])  
                                                        
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
    Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])

    datSpecies=building_tab_pca(InputData_profil,espprinc)

  }

  if (length(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
  }              
  
  
  if(methspecies=="total"){

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
  if (length(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
  }


  if(methspecies=="logevent"){

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
  if (length(objects()[-which(objects() %in% c('methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
  }



# STEP 2 : PCA OR NOT

print("######## STEP 2 PCA/NO PCA ON CATCH PROFILES ########")

  if(acpyesno=="acp"){
    # PCA (Principal Component Analysis)
    marees.pca <- PCA(datSpecies, graph=T)
    
    savePlot(filename='Species projection on the two first factorial axis', type='png', restoreConsole = TRUE)
    dev.off()
    savePlot(filename='Individuals projection on the two first factorial axis', type='png', restoreConsole = TRUE)
    dev.off()
    
    Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   

    graphics.off()

    # Eigenvalues and relative graphics
    marees.pca$eig                                                  

    png("Eigen values.png", width = 1200, height = 800)
    x=1:length(marees.pca$eig[,1])
    barplot(marees.pca$eig[,1],names.arg=x, main="Eigen values")
    dev.off()
    png("Percentage of Inertia.png", width = 1200, height = 800)
    barplot(marees.pca$eig[,2],names.arg=x, main="Percentage of Inertia of factorial axis", xlab="Axis", ylab="% of Inertia")
    dev.off()
    png("Cumulative Percentage of Inertia.png", width = 1200, height = 800) 
    color=rep("grey",length(marees.pca$eig[,1]))
    numaxe=min(which(marees.pca$eig[,3]>70))
    color[1:numaxe]="green"
    barplot(marees.pca$eig[,3],names.arg=x, col=color, main="Cumulative Percentage of Inertia of factorial axis", xlab="Axis", ylab="% of Inertia")
    abline(h=70, col="red")
    text(1,72, "70% of Inertia", col = "red", adj = c(0, -.1))
    dev.off()   
    
    # Data frame given eigenvalues, inertia and cumulative inertia of factorial axis
    tabinertia=data.frame(cbind(Axis=1:length(marees.pca$eig[,1]), Eigenvalues=marees.pca$eig[,1], Inertia=marees.pca$eig[,2], CumulativeInertia=marees.pca$eig[,3]))                  

    # Projection of variables Species on the first factorial axis
    png("Projection of Species on first factorial axis.png", width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(marees.pca,choix="var",axes = c(1, 2),new.plot=FALSE,lim.cos2.var = 0.5)
    plot(marees.pca,choix="var",axes = c(2, 3),new.plot=FALSE,lim.cos2.var = 0.5)
    plot(marees.pca,choix="var",axes = c(1, 3),new.plot=FALSE,lim.cos2.var = 0.5)
    plot(marees.pca,choix="var",axes = c(1, 4),new.plot=FALSE,lim.cos2.var = 0.5)
    plot(marees.pca,choix="var",axes = c(2, 4),new.plot=FALSE,lim.cos2.var = 0.5)
    plot(marees.pca,choix="var",axes = c(3, 4),new.plot=FALSE,lim.cos2.var = 0.5)
    par(op)
    dev.off() 


    # Determine the number of axis
    nbaxes=which(marees.pca$eig[,3]>70)[1]   # we are taking the axis until having 70% of total inertia
    marees.coa=PCA(datSpecies, graph=F, ncp=nbaxes)
    options(digits=6)       # significant digits

    # marees.coa = results of PCA limited to the nbaxes first factorial axis
    datLog=marees.coa$ind$coord
  }

  if (length(objects()[-which(objects() %in% c('methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
  } 



  if(acpyesno=="noacp"){
    datLog=datSpecies
    tabinertia="No PCA"
  }

  if (length(objects()[-which(objects() %in% c('methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
  }




# STEP 3 : CLUSTERING OF THE LOGEVENTS

print("######## STEP 3 CLUSTERING ########")

  if(methmetier=="cah"){
    # Distance between individuals
    marees.dist <- dist(datLog, method=param3)             
    # HAC (Hierarchical Ascendant Classification)
    marees.cah <- hclust(marees.dist, method=param4)

    if (length(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
    }

    # Determine the number of cluster thanks to the scree-test
    inerties.vector=marees.cah$height[order(marees.cah$height,decreasing=T)]
    nb.finalclusters=which(scree(inerties.vector)$epsilon<0)[1]
    nb.finalclusters

    # Cut the dendogram at the selected level
    cah_cluster=cutree(marees.cah,k=nb.finalclusters)

    # Calculation of each cluster size
    n=dim(dat)[1]
    effgroupe=numeric()
    for(k in 1:nb.finalclusters){
      effgroupe[k]=length(which(cah_cluster==k))
    }   
                                                
    # Rectangles plotting
    png("CAHdendogram.png", width = 1200, height = 800)
    plot(marees.cah)
    rect.hclust(marees.cah, k=nb.finalclusters)
    dev.off()

    resval=test.values(cah_cluster,datSpecies)
    target=targetspecies(resval)

    return(list(cah_cluster=cah_cluster, effgroupe=effgroupe, datSpecies=datSpecies, tabinertia=tabinertia, datLog=datLog, nametarget=target$tabnomespcib, namenotarget=target$tabnomesppascib))

  }

  if (length(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
  }



  if(methmetier=="kmeans"){
    # Calculation of optimal k thanks to the variance intraclass
    varintra <- numeric()
    for (k in 2:10){
      clkmeans<-kmeans(datLog, k, iter.max=20, nstart=5)
      varintra [k] <- sum(clkmeans$withinss)
    }
    png("varintrakmeans.png", width = 1200, height = 800)
    plot(varintra)
    dev.off()

    if (length(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
    }

    diffvarintra=diff(varintra,na.rm=T)
    diffdiffvar=diff(diffvarintra,na.rm=T)
    maxdiffdiff=max(diffdiffvar,na.rm=T)
    k=which(diffdiffvar==maxdiffdiff)+1

    # KMEANS
    clkmeans<-kmeans(datLog, k, iter.max=20, nstart=6, algorithm="Hartigan-Wong")

    resval=test.values(clkmeans$cluster,datSpecies)
    target=targetspecies(resval)
    
    #### Inertie

    # Mean profils
    p=length(dat[1,])
    nbclust=length(clkmeans$size)
    mprofil=numeric()
    for(i in 1:nbclust){
      mprofilclusti=mean(dat[which(clkmeans$cluster==i),2:p])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png("Mean profil for each cluster.png", width = 1200, height = 800)
    op <- par(mfrow=c(rep(ceiling(sqrt(nbclust)),2)))
    for(i in 1:nbclust){
      barplot(mprofil[i,])
    }
    par(op)
    title(main="Mean profil for each cluster", xlab="Species", ylab="Euros")
    dev.off()
    
    # Number of Logevents by cluster
    x=c(1:nbclust)
    png("Number of Logevents by cluster.png", width = 1200, height = 800)
    barplot(clkmeans$size, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    dev.off()
    
    # Target Species profils (test-value)
    targetprofil=numeric()
    for(i in 1:nbclust){
      numtargeti=as.numeric(target$tabnumespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      nomtargeti=as.character(target$tabnomespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      targetresvalclusti=resval[numtargeti,i]
      targetresvalclustiord=targetresvalclusti[order(numtargeti,decreasing=F)]
      targetprofili=rep(0,p-1)
      targetprofili[numtargeti]=targetresvalclustiord
      targetprofil=rbind(targetprofil,targetprofili)
    }
    colnames(targetprofil)=colnames(dat[1,2:p])
    png("Target species profil for each cluster.png", width = 1200, height = 800)
    op <- par(mfrow=c(rep(ceiling(sqrt(nbclust)),2)))
    for(i in 1:nbclust){
      barplot(targetprofil[i,])
    }
    par(op)
    title(main="Target species profil for each cluster", xlab="Species", ylab="Test-value")
    dev.off()
    
    return(list(clkmeans=clkmeans, datSpecies=datSpecies, tabinertia=tabinertia, datLog=datLog, nametarget=target$tabnomespcib, namenotarget=target$tabnomesppascib))

  }                                                                                                         

  if (length(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
  }



  if(methmetier=="pam"){

    # Calculation of optimal k thanks to the silhouette
    logpam.silcoeff <- numeric()
    for (k in 2:10){
      logpam=pam(datLog,k)
      logpam.silcoeff [k] <- logpam$silinfo$avg.width
    }

    png("silcoeffpam.png", width = 1200, height = 800)
    plot(logpam.silcoeff)     # k optimal corresponds to maximum of silhouette's coefficients
    dev.off()

    if (length(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
    }

    logpam.silcoeff
    max=max(logpam.silcoeff, na.rm=T)
    k=which(logpam.silcoeff==max)
    logpam=pam(datLog,k)   # PAM with optimal k
    summary(logpam)

#    # Silhouette of PAM with optimal k
#    logpam.si=silhouette(logpam)
#
#    png("silhouetpamopt.png", width = 1200, height = 800)
#    plot(logpam.si)
#    dev.off()

    resval=test.values(logpam$clustering,datSpecies)
    target=targetspecies(resval)
    
    # Mean profils
    p=length(dat[1,])
    nbclust=length(logpam$id.med)
    mprofil=numeric()
    for(i in 1:nbclust){
      mprofilclusti=mean(dat[which(logpam$clustering==i),2:p])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png("Mean profil for each cluster.png", width = 1200, height = 800)
    op <- par(mfrow=c(rep(ceiling(sqrt(nbclust)),2)))
    for(i in 1:nbclust){
      barplot(mprofil[i,])
    }
    par(op)
    title(main="Mean profil for each cluster", xlab="Species", ylab="Euros")
    dev.off()
    
    # Number of Logevents by cluster
    x=c(1:nbclust)
    png("Number of Logevents by cluster.png", width = 1200, height = 800)
    barplot(logpam$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    dev.off()
    
    # Target Species profils (test-value)
    targetprofil=numeric()
    for(i in 1:nbclust){
      numtargeti=as.numeric(target$tabnumespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      nomtargeti=as.character(target$tabnomespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      targetresvalclusti=resval[numtargeti,i]
      targetresvalclustiord=targetresvalclusti[order(numtargeti,decreasing=F)]
      targetprofili=rep(0,p-1)
      targetprofili[numtargeti]=targetresvalclustiord
      targetprofil=rbind(targetprofil,targetprofili)
    }
    colnames(targetprofil)=colnames(dat[1,2:p])
    png("Target species profil for each cluster.png", width = 1200, height = 800)
    op <- par(mfrow=c(rep(ceiling(sqrt(nbclust)),2)))
    for(i in 1:nbclust){
      barplot(targetprofil[i,])
    }
    par(op)
    title(main="Target species profil for each cluster", xlab="Species", ylab="Test-value")
    dev.off()

    return(list(logpam=logpam, datSpecies=datSpecies, tabinertia=tabinertia, datLog=datLog, nametarget=target$tabnomespcib, namenotarget=target$tabnomesppascib))

  }

  if (length(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
  }



  if(methmetier=="clara"){
    n=dim(datLog)[1]

    # Calculation of optimal k thanks to the silhouette
    logclara.silcoeff <- numeric()
    for (k in 2:10){
      logclara=clara(datLog, k, metric=param3, stand=F, samples=5, sampsize=min(n,40+2*k))
      logclara.silcoeff [[k]] <- logclara$silinfo$avg.width
    }

    png("silcoeffclara.png", width = 1200, height = 800)
    plot(logclara.silcoeff)               # k optimal corresponds to maximum of silhouette's coefficients
    dev.off()
    
    if (length(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])   
    }

    logclara.silcoeff
    max=max(logclara.silcoeff, na.rm=T)
    k=which(logclara.silcoeff==max)
    logclara=clara(datLog, k, metric=param3, stand=F, samples=5, sampsize=min(n,40+2*k))  # CLARA with optimal k
    summary(logclara)

#    #Silhouette of CLARA with optimal k
#    logclara.si=silhouette(logclara, full=T)        # full silhouette obtained from the best sample
#
#    png("silhouetclaraopt.png", width = 1200, height = 800)
#    plot(silhouette(logclara, full=T))
#    dev.off()

# Graph cluster
#  clusplot(logclara,col.p = logclara$clustering)

    resval=test.values(logclara$cluster,datSpecies)
    target=targetspecies(resval)

    # Mean profils
    nbSpec=length(datSpecies[1,])
    nbclust=length(logclara$i.med)
    mprofil=numeric()
    for(i in 1:nbclust){
      mprofilclusti=mean(datSpecies[which(logclara$clustering==i),])
      mprofil=rbind(mprofil,mprofilclusti)
    }
    png("Mean profil for each cluster.png", width = 1200, height = 800)
    op <- par(mfrow=c(rep(ceiling(sqrt(nbclust)),2)))
    for(i in 1:nbclust){
      barplot(mprofil[i,])
    }
    par(op)
    title(main="Mean profil for each cluster", xlab="Species", ylab="Euros")
    dev.off()
    
    # Number of Logevents by cluster
    x=c(1:nbclust)
    png("Number of Logevents by cluster.png", width = 1200, height = 800)
    barplot(logclara$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    dev.off()
    
    # Target Species profils (test-value)
    targetprofil=numeric()
    targetresval=numeric()
    for(i in 1:nbclust){
      numtargeti=as.numeric(target$tabnumespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      nomtargeti=as.character(target$tabnomespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      targetresvalclusti=rep(0,nbSpec)
      targetresvalclusti[numtargeti]=resval[numtargeti,i]
      #targetresvalclusti=resval[numtargeti,i]
      targetresval=rbind(targetresval,targetresvalclusti)
#      targetresvalclustiord=targetresvalclusti[order(numtargeti,decreasing=F)]
#      targetprofili=rep(0,nbSpec)
#      targetprofili[numtargeti]=targetresvalclustiord
#      targetprofil=rbind(targetprofil,targetprofili)
    }
    #colnames(targetprofil)=colnames(datSpecies)
    colnames(targetresval)=colnames(datSpecies)
    png("Target species profil for each cluster.png", width = 1200, height = 800)
    op <- par(mfrow=c(rep(ceiling(sqrt(nbclust)),2)))
    for(i in 1:nbclust){
      #barplot(targetprofil[i,])
      barplot(targetresval[i,])
    }
    par(op)
    title(main="Target species profil for each cluster", xlab="Species", ylab="Test-value")
    dev.off()
    
    return(list(logclara=logclara, datSpecies=datSpecies, tabinertia=tabinertia, datLog=datLog, nametarget=target$tabnomespcib, namenotarget=target$tabnomesppascib))

  }
  
  if (length(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])>0){
      Store(objects()[-which(objects() %in% c('dat','methspecies','param1','param2','acpyesno','methmetier','param3','param4'))])  
  }
  # end of the methods

} # end of the function "classif"

