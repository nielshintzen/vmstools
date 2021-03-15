################################################################################
#  STEP 3 OF THE MULTIVARIATE CLASSIFICATION :                                 #
#         RUN THE CLUSTERING OF THE LOGEVENTS                                  #
#         4 METHODS ARE AVALAIBLE : HAC / KMEANS / PAM / CLARA                 #
################################################################################




#' Finding metiers from a reduced EFLALO dataset, step 3: clustering logevents
#' using various multivariate methods
#' 
#' This function represents the third step in the workflow processing logbooks
#' data for identifying metiers.
#' 
#' This step allows applying various clustering analyses on the data sets
#' coming out of the first and second step. All methods will lead to a
#' classification of all individuals (logevents), but they differ in their
#' nature and then consequently in their outcomes. The four methods available
#' are - Hierarchical Ascending Classification (HAC), with user-defined method
#' for estimating distances and link for aggregating individuals\cr -
#' K-Means,\cr - Partitioning Around Medoids (PAM),\cr - Clustering LARge
#' Applications (CLARA).
#' 
#' The HAC method works by calculating the distance between individuals using
#' the method selected with param1 ("euclidean", "maximum", "manhattan",
#' "canberra", "binary" or "minkowski") and aggregating them based on the
#' distance between clusters, using the link selected with param2 ("ward",
#' "single", "complete", "average", "mcquitty", "median" or "centroid"). In
#' HAC, the number of classes is determined afterwards, once all combinations
#' have been calculated, by using the objective criteria of scree test which
#' detects the third marginal loss of inertia between two consecutive numbers
#' of classes. Therefore, the computing time and memory request for this method
#' can be quite comprehensive, and may reach memory limits on standard PC when
#' operating with very large datasets.
#' 
#' The K-Means method works by randomly choosing k individuals, or kernels (k
#' corresponding to the final number of classes), and then affecting each
#' individuals of the dataset to the closest kernel. Each time, the gravity
#' center of the class is recalculated, thus reinitialising the calculation of
#' distances to the next individual. In order to define the most appropriate
#' number of classes, this procedure is repeated with differents values for k,
#' from 2 to 15. The final number of classes is identified by using the
#' criteria of scree test which detects the third marginal loss of inertia
#' between two consecutive numbers of classes.
#' 
#' The PAM method works slightly around the same principle, starting with the
#' initialisation of k medoids. The medoid is the individual in a class which
#' shows least dissimilarity with other individuals in the same class, and the
#' remaining individuals are affected to their closest medoid. Then the sum of
#' dissimilarities is calculated and compared with the sum of dissimilarities
#' if any other individual in the class had been playing the role of the
#' medoid, and then the medoid is eventually adjusted accordingly, until full
#' stabilisation of the procedure. The most appropriate number of classes is
#' identified by using the estimated silhouette of the classification for each
#' value of k. The silhouette represent an average comparison of the distance
#' between an individual and the other individuals from its class, and between
#' the same individual and the other individuals from the next closest class,
#' and is therefore an objective measurement of the quality of the
#' classification. The final number of classes retained is the one for which
#' the second maximum of the silhouettes is reached. It is to be noted that the
#' PAM method is not designed for working with very large datasets (>10 000
#' lines), and may quickly reach memory limits.
#' 
#' The CLARA method is an extension of the PAM algorithm aiming at working with
#' large datasets, and is therefore more efficient and powerful than PAM.  It
#' works by sampling 5 subsets of the dataset for each value of k and running
#' the PAM algorithm explained above on these subsets, and then to keep only
#' the subset giving the best classification (silhouette). As in PAM, the final
#' number of classes retained is the one for which the second maximum of the
#' silhouettes is reached. Afterwards, all remaining individuals are affected
#' to their closest medoid, using user-defined method (param1) for calculating
#' distances ('euclidean' being used as the default, but 'manhattan' could also
#' be used).
#' 
#' 
#' @param datSpecies numerical matrix with logevents as rows and species as
#' columns, with percentage values (between 0 and 100) of each species in the
#' logevent catches.  Logevent ID (LE_ID) should be as row names. Typically,
#' this input table will be produced from the step 1 of the metier analysis
#' applied on the eflalo initial data, using the function
#' extractTableMainSpecies()
#' @param datLog numerical matrix with logevents as rows, and values to be used
#' for calculating distances between individuals as columns.  Typically, this
#' input table is produced by the step 2 of the metier analysis, using the
#' function getTableAfterPCA(). If a PCA was run, selected Principal Components
#' will appear as columns. If no PCA was run, the matrix will be the same as
#' datSpecies, with percentage values by species.
#' @param analysisName character, the name of the run. Used for the file name
#' of the plots.
#' @param methMetier character. The name of the clustering method to be used.
#' Must be chosen between "hac", "kmeans", "pam" and "clara".
#' @param param1 character. Parameter used for chosing the method calculating
#' distances between individuals, to be used in HAC and CLARA algorithms.  For
#' HAC, it can be chosen between "euclidean", "maximum", "manhattan",
#' "canberra", "binary" or "minkowski". For CLARA, between "euclidean" and
#' "manhattan".  For PAM and K-means, this must be set to NULL.
#' @param param2 character. Parameter used for chosing the method calculating
#' distances between clusters, to be used in HAC algorithm.  For HAC, it can be
#' chosen between "ward", "single", "complete", "average", "mcquitty", "median"
#' or "centroid".  For PAM, CLARA and K-means, this must be set to NULL.
#' @return For any of the four methods used, a number of graphs are produced
#' and saved directly in the working directory.  They describe 1) projections
#' of results on factorial plans, 2) the mean profile of the whole dataset,
#' i.e. the average percentage of the various species in a logevent across all
#' individuals, 3) the mean and standard deviation profile in terms of average
#' percentage of catch by species within each cluster, 4) the number of
#' logevents by clusters, and 5) the profile of test-values by cluster.
#' 
#' Finally, the function returns a list with a number of results and
#' diagnostics on the performance of the method: \item{LE_ID_clust}{a data
#' frame with two columns, linking the initial ID name of the Logevent (LE_ID)
#' with the cluster metier where the ID has been allocated. }
#' \item{clusters}{diagnostics of the clustering process. It may vary between
#' the four methods.  } \item{betweenVarClassifOnTot}{percentage of variance
#' explained by the classification.  } \item{nbClust}{final number of clusters
#' retained.  } \item{summaryClusters}{array documenting, for each cluster, the
#' minimum, mean, maximum, as well as the 25\%, 50\% and 75\% quantiles values
#' of the percentage of catch by species for the individual logevents in the
#' cluster.  } \item{testValues}{matrix of test-values by species and cluster.
#' The test-value measures for each species the difference between the average
#' percentage of catches in the cluster compared to the average percentage of
#' catch in the total dataset, thus large positive values (>1.98) will point
#' out the most characteristic species in the clusters.  }
#' \item{testValuesSpecies}{a tabulated list ranking the most characteristic
#' species in the clusters (ordering species with a test-value > 1.98 by
#' decreasing order of test-value).  } \item{percLogevents}{a matrix giving the
#' percentage of logevents catching each species by cluster. }
#' \item{descClusters}{a data frame giving some descriptive statistics for each
#' cluster, like cluster size, number of species needed to have at least 50\%
#' of the cluster's total catch, number of species with a test-value > 1.98 in
#' the cluster, number of species caught in at least 50\% of the logevents,
#' etc... } \item{tabClusters}{a 3d-array giving a table summurazing for each
#' cluster the most important species (in terms of catch), the associated
#' test-value, and the percentage of logevents of the cluster catching these
#' species.  } \item{targetSpecies}{a list giving the target species by
#' cluster. }
#' @note A number of libraries are initially called for the whole metier
#' analyses and must be installed :
#' (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)
#' @author Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER,
#' France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)
#' @seealso \code{\link{getEflaloMetierLevel7}},
#' \code{\link{selectMainSpecies}}, \code{\link{extractTableMainSpecies}},
#' \code{\link{getMetierClusters}}, \code{\link{getTableAfterPCA}}
#' @references Development of tools for logbook and VMS data analysis. Studies
#' for carrying out the common fisheries policy No MARE/2008/10 Lot 2
#' @examples
#' 
#' 
#'   \dontrun{
#'  
#'   data(eflalo)
#'   
#'   eflalo <- formatEflalo(eflalo)
#' 
#'   eflalo <- eflalo[eflalo$LE_GEAR=="OTB",]
#'   # note that output plots will be sent to getwd()
#'   analysisName <- "metier_analysis_OTB"
#'   
#'   dat <- eflalo[,c("LE_ID",grep("EURO",colnames(eflalo),value=TRUE))]
#'   names(dat)[-1] <- unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))
#' 
#'   explo <- selectMainSpecies(dat, analysisName, RunHAC=TRUE, DiagFlag=FALSE)
#'     #=> send the LE_ID and LE_KG_SP columns only
#'            
#'   Step1 <- extractTableMainSpecies(dat, explo$NamesMainSpeciesHAC,
#'               paramTotal=95, paramLogevent=100)
#'     #=> send the LE_ID and LE_KG_SP columns only             
#' 
#'   rowNamesSave <- row.names(Step1)
#'   row.names(Step1) <- 1:nrow(Step1)
#' 
#'   # Run a PCA
#'   Step2 <- getTableAfterPCA(Step1, analysisName, pcaYesNo="pca",
#'               criterion="70percents")
#' 
#'   row.names(Step1) <- rowNamesSave
#'   row.names(Step2) <- rowNamesSave
#'   
#'   # Define a metier for each logevent running the CLARA algorithm 
#'   Step3 <- getMetierClusters(Step1, Step2, analysisName,
#'               methMetier="clara", param1="euclidean", param2=NULL)
#'               
#'   }
#' 
#' 
#' @export getMetierClusters
getMetierClusters = function(datSpecies,datLog,analysisName="",methMetier="clara",param1="euclidean",param2=NULL){

  # Load the table linking 3A-CODE (FAO CODE of species) to the species assemblage (level 5).
  data(correspLevel7to5)
  require(lattice)

  LE_ID=rownames(datSpecies)
  nbSpec=ncol(datSpecies)
  datSpecies=as.matrix(datSpecies,ncol=nbSpec,nrow=length(LE_ID))

  print("######## STEP 3 CLUSTERING ########")

  t1=Sys.time()
  print(paste(" --- selected method :",methMetier, " ---"))


########################################################################################################################################   HAC

  if(methMetier=="hac"){

    classifWithinVar=numeric()
    classifBetweenVar=numeric()
    classifQuality=numeric()
    sampleList=numeric()
    mProfilSample=numeric()
    classifVarExplain=numeric()

    nbLog=nrow(datLog)
    nbDim=ncol(datLog)

    # Center of gravity of datLog
    centerOfGravityDatLog=colMeans(datLog)

    # HAC like CLARA (HAC on sample, affectation of each logevent to a cluster, quality of classification, do it 5 times, choose the sample which gives the best quality of classification)
    print("hac on subsets...")

    for(i in 1:5){

      numSample=i
      print(paste("sample",i))
      # Sample of size 10000 logevents or 30% of all logevents
      minsam=min(nbLog,max(10000,round(nbLog*30/100)))
      sam=sample(1:nbLog,size=minsam,replace=FALSE)
      # Record the 5 samples
      sampleList=rbind(sampleList,sam)
      outofsam=setdiff(1:nbLog,sam)
      sampleDatLog=datLog[sam,]
      sampleDatSpecies=datSpecies[sam,]

      # HAC on the sample
      log.hac=hcluster(sampleDatLog, method=param1, link=param2)
      inerties.vector=log.hac$height[order(log.hac$height,decreasing=TRUE)]
      nbClust=which(scree(inerties.vector)[,"epsilon"]<0)[3]

      # Cut the dendogram at the selected level
      sampleClusters=cutree(log.hac,k=nbClust)

      Store(objects())
      gc(reset=TRUE)

      # Add the cluster to each logevent of the sample
      sampleDatLogWithClusters=cbind(sampleDatLog,sampleClusters)

      sampleClusters=sampleDatLogWithClusters[,ncol(sampleDatLogWithClusters)]

      # Within and between variance of clusters and classification
      centerOfGravityClassif=numeric()
      withinVarClusters=numeric()
      sizeClusti=numeric()
      centerOfGravitySampleDatLog=colMeans(sampleDatLog)
      centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravitySampleDatLog)
      for(k in 1:nbClust){  # Within variance by cluster

        clusti=sampleDatLogWithClusters[which(sampleClusters==k),1:nbDim]
        if(length(which(sampleClusters==k))==1)  centerOfGravityClusti=clusti
        else centerOfGravityClusti=colMeans(clusti)
        centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
        sizeClusti[k]=length(which(sampleClusters==k))
        if(length(which(sampleClusters==k))==1)  withinVarClusters[k]=0
        else withinVarClusters[k]=sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))

      }
      # Between variance
      classifBetweenVar=(1/nbLog)*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
      # Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
      withinVarClusterOnTot=(1/nbLog)*sum(withinVarClusters)/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
      betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
      classifVarExplain=c(classifVarExplain,betweenVarClassifOnTot)


      # Catch profile by cluster for each sample
      nbSpec=ncol(datSpecies)
      mprofil=numeric()
      blank=rep(00,nbSpec)
      for(k in 1:nbClust){
        mprofilclusti=colMeans(sampleDatSpecies[which(sampleClusters==k),])
        mprofil=rbind(mprofil,mprofilclusti)
      }
      mprofil=rbind(mprofil,blank)

      mProfilSample=rbind(mProfilSample,mprofil)


      # Graphics

      # Calculation of each cluster size
      sizeClusters=numeric()
      for(k in 1:nbClust){
        sizeClusters[k]=length(which(sampleClusters==k))
      }

      # Compute the test-values for species
      resval=test.values(sampleClusters,sampleDatSpecies)
      # Determine the target species
      target=targetspecies(resval)


      # Projections on the first factorial plans
      png(paste(analysisName,numSample,"Sample_Projections.png",sep="_"), width = 1200, height = 800)
      op <- par(mfrow=c(2,3))
      plot(sampleDatLog[,1], sampleDatLog[,2], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
      if(dim(datLog)[2]>2) {
      plot(sampleDatLog[,2], sampleDatLog[,3], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
      plot(sampleDatLog[,1], sampleDatLog[,3], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
      if(dim(datLog)[2]>3) {
      plot(sampleDatLog[,1], sampleDatLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
      plot(sampleDatLog[,2], sampleDatLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
      plot(sampleDatLog[,3], sampleDatLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(sampleClusters)], main="Projection of HAC classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
      }}
      par(op)
      dev.off()


      # Rectangles plotting
      png(paste(analysisName,numSample,"Sample_Dendogram.png",sep="_"), width = 1200, height = 800)
      plclust(log.hac,labels=FALSE,hang=-1,ann=FALSE)
      title(main="HAC dendogram",xlab="Logevents",ylab="Height")
      rect.hclust(log.hac, k=nbClust)
      dev.off()


      # Catch profile of the dataset
      meanprofile=colMeans(sampleDatSpecies)
      png(paste(analysisName,numSample,"Sample_Catch profile of the sample.png",sep="_"), width = 1200, height = 800)
      op <- par(las=2)
      barplot(meanprofile, main="Catch profile of the sample", xlab="Species", ylab="Percentage of catch")
      par(op)
      mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
      dev.off()


      # Catch profile by cluster
      nbSpec=ncol(sampleDatSpecies)
      summarySampleClusters=array(0,dim=c(6,nbSpec,nbClust))
      dimnames(summarySampleClusters)[[1]]=c("Min.","1st Qu.","Median", "Mean", "3rd Qu.", "Max.")
      dimnames(summarySampleClusters)[[2]]=names(meanprofile)
      dimnames(summarySampleClusters)[[3]]=paste("Cluster",1:nbClust)
      for(k in 1:nbClust){
        if(sizeClusters[k]==1){
          summarySampleClusters[,,k]=apply(t(as.matrix(sampleDatSpecies[which(sampleClusters==k),])),2,
            function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
        }else{
          summarySampleClusters[,,k]=apply(sampleDatSpecies[which(sampleClusters==k),],2,
            function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
        }
      }
      # Species names for catch profile plots
      nameSpPlot=character()
      catchMeanThreshold=2
      for(k in 1:nbClust){
        namSpi=names(meanprofile[which(t(summarySampleClusters["Mean",,k])>catchMeanThreshold)])
        numSpi=which(t(summarySampleClusters["Mean",,k])>catchMeanThreshold)
        nameSpPloti=rep("",nbSpec)
        nameSpPloti[numSpi]=namSpi
        nameSpPlot=rbind(nameSpPlot,nameSpPloti)
      }
      png(paste(analysisName,numSample,"Sample_Mean profile by cluster of the sample.png",sep="_"), width = 1200, height = 800)
      op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
      for(k in 1:nbClust){
        op2 <- par(las=2)
        barplot(t(summarySampleClusters["Mean",,k]), names.arg=nameSpPlot[k,], xlab="Species", ylab="Percentage of catch", col="gray")
        par(op2)
        mtext(paste("Cluster",k), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
      }
      par(op)
      title(main=paste("Catch profile by cluster of the sample","\n","\n",sep=""))
      dev.off()


      # Standard deviation profile by cluster
      sdprofil=matrix(0,nrow=nbClust,ncol=nbSpec)
      namSdPlot=character()
      SdThreshold=2
      for(k in 1:nbClust){
        if(length(which(sampleClusters==k))==1){ sdprofilclusti=rep(0,nbSpec)
        }else{sdprofilclusti=apply(sampleDatSpecies[which(sampleClusters==k),],2,sd)}
        namSDi=names(which(sdprofilclusti>SdThreshold))
        numSDi=which(sdprofilclusti>SdThreshold)
        namSdPloti=rep("",nbSpec)
        namSdPloti[numSDi]=namSDi
        sdprofil[k,]=sdprofilclusti
        namSdPlot=rbind(namSdPlot,namSdPloti)
      }
      rownames(sdprofil) <- 1:nrow(sdprofil)
      png(paste(analysisName,numSample,"Sample_Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
      op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
      for(k in 1:nbClust){
        op2 <- par(las=2)
        barplot(sdprofil[k,], names.arg=namSdPlot[k,], xlab="Species", ylab="Percentage of catch")
        par(op2)
        mtext(paste("Cluster",k), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
      }
      par(op)
      title(main=paste("Standard deviation profile by cluster","\n","\n",sep=""))
      dev.off()


      # Number of Logevents by cluster
      x=c(1:nbClust)
      png(paste(analysisName,numSample,"Sample_Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
      coord=barplot(sizeClusters, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
      barplot(sizeClusters, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
      text(coord,sizeClusters+100,sizeClusters,font=2,xpd=NA)
      dev.off()


      # Target Species profiles (test-value)
      targetresval=numeric()
      nameTargetPlot=character()
      for(k in 1:nbClust){
        nomtargeti=as.character(target$tabnomespcib[k,which(!is.na(target$tabnumespcib[k,]))])
        numtargeti=as.numeric(target$tabnumespcib[k,which(!is.na(target$tabnumespcib[k,]))])
        nameTargetPloti=rep("",nbSpec)
        nameTargetPloti[numtargeti]=nomtargeti
        nameTargetPlot=rbind(nameTargetPlot,nameTargetPloti)
        targetresvalclusti=rep(0,nbSpec)
        targetresvalclusti[numtargeti]=resval[nomtargeti,k]
        targetresval=rbind(targetresval,targetresvalclusti)
      }

      png(paste(analysisName,numSample,"Sample_Profile of target species by cluster.png",sep="_"), width = 1200, height = 800)
      op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
      for(k in 1:nbClust){
        op2 <- par(las=2)
        barplot(targetresval[k,],names.arg=nameTargetPlot[k,], cex.names=1, xlab="Species", ylab="Test-value")
        par(op2)
        mtext(paste("Cluster",k), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
      }
      par(op)
      title(main=paste("Profile of target species by cluster","\n","\n",sep=""))
      dev.off()

      Store(objects())
      gc(reset=TRUE)

    } # end of for(i in 1:5)




    # Select the sample which gives the smaller classification's quality (the best sample)
    sam=sampleList[which.max(classifVarExplain),]
    outofsam=setdiff(1:nbLog,sam)
    sampleDatLog=datLog[sam,]

    nbLogSample=nrow(sampleDatLog)
    nbDim=ncol(sampleDatLog)


    # HAC with the best sample
    print("Final HAC")
    log.hac=hcluster(sampleDatLog, method=param1, link=param2)


    # Determine the number of cluster thanks to the scree-test
    inerties.vector=log.hac$height[order(log.hac$height,decreasing=TRUE)]
    nbClust=which(scree(inerties.vector)[,"epsilon"]<0)[3]

    # Cut the dendogram at the selected level
    sampleClusters=cutree(log.hac,k=nbClust)
    sampleClusters=as.factor(sampleClusters)

    sampleDatLogWithClusters=data.frame()
    sampleDatLogWithClusters=cbind(sampleDatLog,sampleClusters)
    sampleDatLogWithClusters=as.data.frame(sampleDatLogWithClusters)

    # Discriminante analysis on the learning dataset
    learning=lda(sampleClusters~.,data=sampleDatLogWithClusters)

    otherLog=datLog[outofsam,]
    otherLog=as.data.frame(otherLog)

    # Predict the cluster for the other logevent
    pred=predict(learning,otherLog)
    otherDatLogWithClusters=cbind(otherLog,pred$class)
    colnames(otherDatLogWithClusters)=colnames(sampleDatLogWithClusters)

    # Rebuilt complete datLog with clusters
    clusters=numeric(length=nbLog)
    clusters[sam]=sampleClusters
    clusters[outofsam]=pred$class


    # Within and between variance of clusters and classification
    centerOfGravityClassif=numeric()
    withinVarClusters=numeric()
    sizeClusti=numeric()
    centerOfGravityDatLog=colMeans(datLog)
    centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)
    for(k in 1:nbClust){  # Within variance by cluster

      clusti=datLog[which(clusters==k),1:nbDim]
      if(length(which(clusters==k))==1)  centerOfGravityClusti=clusti
      else centerOfGravityClusti=colMeans(clusti)
      centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
      sizeClusti[k]=length(which(clusters==k))
      if(length(which(clusters==k))==1)  withinVarClusters[k]=0
      else withinVarClusters[k]=sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))

    }
    # Between variance
    classifBetweenVar=(1/nbLog)*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
    # Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
    withinVarClusterOnTot=(1/nbLog)*sum(withinVarClusters)/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
    betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100


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
    rownames(target$tabnomespcib)=paste("Cluster",1:nbClust)


    # Compute the percentage of logevents catching each species by cluster
    mainSpecies=colnames(datSpecies)
    percLogevents=matrix(0,ncol=length(mainSpecies),nrow=nbClust,dimnames=list(paste("Cluster ",1:nbClust,sep=""),mainSpecies))
    for(i in 1:nbClust){
      percLogevents[i,]=round(sapply(mainSpecies,function(x) (sizeClusters[i]-length(which(datSpecies[clusters==i,x]==0)))/sizeClusters[i]*100),digits=1)
    }
    

    # Projections on the first factorial plans
    png(paste(analysisName,"Projections.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
    if(dim(datLog)[2]>2) {
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
    plot(datLog[,1], datLog[,3], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
    if(dim(datLog)[2]>3) {
    plot(datLog[,1], datLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
    plot(datLog[,2], datLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
    plot(datLog[,3], datLog[,4], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="Projection of HAC classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
    }}
    par(op)
    dev.off()

    # For a paper
    X11(5,5)
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="", xlab="Axis 1", ylab="Axis 2")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_1_2_HAC',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="", xlab="Axis 1", ylab="Axis 2")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_1_2_HAC',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(sizeClusters))[as.numeric(clusters)], main="", xlab="Axis 2", ylab="Axis 3")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_2_3_HAC',sep="_"), type='png')
    dev.off()


    # Rectangles plotting
    png(paste(analysisName,"Dendogram.png",sep="_"), width = 1200, height = 800)
    plclust(log.hac,labels=FALSE,hang=-1,ann=FALSE)
    title(main="HAC dendogram",xlab="Logevents",ylab="Height")
    rect.hclust(log.hac, k=nbClust)
    dev.off()


    # Catch profile of the dataset
    meanprofile=colMeans(datSpecies)
    png(paste(analysisName,"Catch profile of the dataset.png",sep="_"), width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Catch profile of the dataset", xlab="Species", ylab="Percentage of catch")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    dev.off()


    # Catch profile by cluster
    nbSpec=ncol(datSpecies)
    summaryClusters=array(0,dim=c(6,nbSpec,nbClust))
    dimnames(summaryClusters)[[1]]=c("Min.","1st Qu.","Median", "Mean", "3rd Qu.", "Max.")
    dimnames(summaryClusters)[[2]]=names(meanprofile)
    dimnames(summaryClusters)[[3]]=paste("Cluster",1:nbClust)
    for(i in 1:nbClust){
      if(sizeClusters[i]==1){
        summaryClusters[,,i]=apply(t(as.matrix(datSpecies[which(clusters==i),])),2,
          function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
      }else{
        summaryClusters[,,i]=apply(datSpecies[which(clusters==i),],2,
          function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
      }
    }
    # Species names for catch profile plots
    nameSpPlot=character()
    catchMeanThreshold=2
    for(i in 1:nbClust){
      namSpi=names(meanprofile[which(t(summaryClusters["Mean",,i])>catchMeanThreshold)])
      numSpi=which(t(summaryClusters["Mean",,i])>catchMeanThreshold)
      nameSpPloti=rep("",nbSpec)
      nameSpPloti[numSpi]=namSpi
      nameSpPlot=rbind(nameSpPlot,nameSpPloti)
    }
    png(paste(analysisName,"Catch profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(t(summaryClusters["Mean",,i]), names.arg=nameSpPlot[i,], xlab="Species", ylab="Percentage of catch", col="gray")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Catch profile by cluster","\n","\n",sep=""))
    dev.off()


    # For a paper : levelplot
    X11(4,4)
    mat <- t(summaryClusters["Mean",,])
    #rownames(mat) <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
    #                          "XVI","XVII","XVIII","XIX","XX")[1:nrow(mat)]
    rownames(mat) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")[1:nrow(mat)]
    sp <- apply(mat,2, function(x) length(which(x>20)))
    sp2=sp[which(sp>=1)]
    #colnames(mat)[sp==0] <- ""
    mat2=mat[,which(colnames(mat)%in%names(sp2))]
    cc <- colorRampPalette(c("white", "black"),space = "rgb", interpolate="spline")
    print(levelplot(mat2, cuts=4, aspect=1, xlab="", ylab="", col.regions=cc(5), at=c(0,20,40,60,80,100), scales=list(cex=0.7), colorkey=list(space="right", at=c(0,20,40,60,80,100), width=1.1)))
    savePlot(filename=paste(analysisName,'mean_profile_by_cluster_levelplot',sep="_"), type='png')
    dev.off()


    # Standard deviation profile by cluster
    sdprofil=matrix(0,nrow=nbClust,ncol=nbSpec)
    namSdPlot=character()
    SdThreshold=2
    for(i in 1:nbClust){
      if(length(which(clusters==i))==1){ sdprofilclusti=rep(0,nbSpec)
      }else{sdprofilclusti=apply(datSpecies[which(clusters==i),],2,sd)}
      namSDi=names(which(sdprofilclusti>SdThreshold))
      numSDi=which(sdprofilclusti>SdThreshold)
      namSdPloti=rep("",nbSpec)
      namSdPloti[numSDi]=namSDi
      sdprofil[i,]=sdprofilclusti
      namSdPlot=rbind(namSdPlot,namSdPloti)
    }
    rownames(sdprofil) <- 1:nrow(sdprofil)
    png(paste(analysisName,"Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(sdprofil[i,], names.arg=namSdPlot[i,], xlab="Species", ylab="Percentage of catch")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Standard deviation profile by cluster","\n","\n",sep=""))
    dev.off()


    # Number of Logevents by cluster
    x=c(1:nbClust)
    png(paste(analysisName,"Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
    coord=barplot(sizeClusters, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    barplot(sizeClusters, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
    text(coord,sizeClusters+400,sizeClusters,font=2,xpd=NA)
    dev.off()


    # Profile of test-values by cluster
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

    png(paste(analysisName,"Profile of test-values by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,],names.arg=nameTargetPlot[i,], cex.names=1, xlab="Species", ylab="Test-value")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Profile of test-values by cluster","\n","\n",sep=""))
    dev.off()




    # Descriptive tables of the clusters
    clusterDesc=matrix(0,nrow=9,ncol=nbClust)
    for(i in 1:nbClust){
      clusterDesc[,i]=c(i,
                        length(which(cumsum(t(summaryClusters["Mean",,i])[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)])<50))+1,
                        length(which(cumsum(t(summaryClusters["Mean",,i])[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)])<90))+1,
                        length(which(t(summaryClusters["Median",,i])>50)),
                        length(which(resval[,i]>1.96)),
                        length(which(resval[,i]>3.29)),
                        length(which(apply(datSpecies,2,function (x) (sizeClusters[i]-length(which(x[clusters==i]==0)))/sizeClusters[i]*100)>50)),
                        length(which(apply(datSpecies,2,function (x) (sizeClusters[i]-length(which(x[clusters==i]==0)))/sizeClusters[i]*100)>90)),
                        sizeClusters[i])
    }
    rownames(clusterDesc)=c("Number of species",
                            "to have 50% of catch", "to have 90% of catch",
                            "with a median higher than 50",
                            "with a test-value > 1.96", "with a test-value > 3.29",
                            "catch in 50% of the logevents", "catch in 90% of the logevents",
                            "Clusters size")
    colnames(clusterDesc)=1:nbClust
    clusterDesc2=as.data.frame(clusterDesc)


    # Summary tables of the clusters
    namesSpecies=matrix(NA,nrow=nbClust,ncol=10)
    namesCapt=matrix(NA,nrow=nbClust,ncol=5)
    nbSpeciesCatch = min(5,dim(t(summaryClusters["Mean",,]))[[2]])
    namesTarget=matrix(NA,nrow=nbClust,ncol=5)
    nbSpeciesVT = min(5,dim(target$tabnomespcib)[[2]])
    tabLibname=matrix(NA,nrow=nbClust,ncol=10)
    listLibname=list()
    
    for(i in 1:nbClust){
      namesCapt[i,]=colnames(t(summaryClusters["Mean",,i]))[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)][1:nbSpeciesCatch]
      a=as.data.frame(t(summaryClusters["Mean",target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])],i]))
      colnames(a)= target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])]
      if(length(a)!=0){
        namesTarget[i,1:length(target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])])]=colnames(a[order(a,decreasing=TRUE)])
      }
      namesSpecies[i,1:length(union(namesCapt[i,],namesTarget[i,]))]=union(namesCapt[i,],namesTarget[i,])
    }
    
    for(i in 1:nbClust){
        listLibname[[i]]=lapply(as.list(namesSpecies[i,]), function(x) if(length(which(correspLevel7to5[,"X3A_CODE"]==x))==0) "NA" 
                                                    else correspLevel7to5[which(correspLevel7to5[,"X3A_CODE"]==x),"French_name"])
        tabLibname[i,]=unlist(lapply(listLibname[[i]], function(x) as.character(unlist(x))))
    }

    tabPropCatch=matrix(NA,nrow=nbClust,ncol=10)
    tabTestVal=matrix(NA,nrow=nbClust,ncol=10)
    tabPropLog=matrix(NA,nrow=nbClust,ncol=10)

    for(i in 1:nbClust){
      print("-----------------------------------------------------------------")
      print(paste("Cluster",i))
      propCatch=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)>=0.1)]
      tabPropCatch[i,1:length(propCatch)]=propCatch
      print(propCatch)
      testVal=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) resval[x,i]),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)>=0.1)]
      tabTestVal[i,1:length(testVal)]=testVal
      print(testVal)
      propLog=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) (sizeClusters[i]-length(which(datSpecies[clusters==i,x]==0)))/sizeClusters[i]*100),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) (sizeClusters[i]-length(which(datSpecies[clusters==i,x]==0)))/sizeClusters[i]*100),digits=1)>=0.1)]
      tabPropLog[i,1:length(propLog)]=propLog
      print(propLog)
    }

    tabClusters=array(0,dim=c(10,5,nbClust))
    dimnames(tabClusters)[[2]]=c("Libname","FAO","Test-value","% Catch","% Logevents")
    dimnames(tabClusters)[[3]]=paste("Cluster",1:nbClust)
    for(i in 1:nbClust){
      tabClusters[,,i]=cbind(tabLibname[i,],namesSpecies[i,],tabTestVal[i,],tabPropCatch[i,],tabPropLog[i,])
    }

    sizeTabClusters=numeric()
    for(i in 1:nbClust){
      sizeTabClusters[i]=min(length(namesSpecies[i,!is.na(namesSpecies[i,])]),length(tabPropCatch[i,!is.na(tabPropCatch[i,])]),length(tabTestVal[i,!is.na(tabTestVal[i,])]),length(tabPropLog[i,!is.na(tabPropLog[i,])]))
    }


    # Target Species
    # Intersection of species from tabClusters having : - % Cumulated Catch > thresholdCatch
    #                                                   - Test-value > thresholdTestValue
    #                                                   - % Logevents > thresholdLogevents 
    thresholdCatch=75
    thresholdTestValue=3
    thresholdLogevents=30
    
    sppCumCatch=list()
    sppTestValue=list()
    sppLogevents=list()
    targetSpeciesByCluster=list()
    
    for (i in 1:nbClust){
      percCatchCum=cumsum(as.numeric(tabClusters[,"% Catch",i]))
      nbSpSel=length(which(percCatchCum<thresholdCatch))+1
      sppCumCatch[[i]]=tabClusters[1:nbSpSel,"FAO",i]
      
      sppTestValue[[i]]=tabClusters[which(as.numeric(tabClusters[,"Test-value",i])>thresholdTestValue),"FAO",i]
      
      sppLogevents[[i]]=tabClusters[which(as.numeric(tabClusters[,"% Logevents",i])>thresholdLogevents),"FAO",i]
    
      targetSpeciesByCluster[[i]]=intersect(sppCumCatch[[i]],sppTestValue[[i]])
      targetSpeciesByCluster[[i]]=intersect(targetSpeciesByCluster[[i]],sppLogevents[[i]])
    }

    # List of metiers (level 7)
    listMetiersL7=list()
    for (i in 1:nbClust){
      metiersClusteri=targetSpeciesByCluster[[i]]
      metiersClusteri=as.character(unique(unlist(metiersClusteri)))
      metiersClusteri=paste(unlist(strsplit(metiersClusteri," ")),collapse=" ")
      listMetiersL7[[i]]=metiersClusteri
    }
  
    # Metier (level 7) of each logevent
    metierByLogeventL7=unlist(sapply(clusters,function(x) listMetiersL7[[x]]))



    # Create csv tables
    write.table(clusterDesc2,file="descClusters.csv",col.names=FALSE,sep=";")
    
    dfClust=data.frame()
    dfClust=paste("Clust ",1:nbClust,sep="")
    for(i in 1:nbClust){
      write.table(dfClust[i],file="tabClusters.csv",append=TRUE,row.names=FALSE,sep=";")
      tabClusti=as.data.frame(tabClusters[1:sizeTabClusters[i],,i])
      write.table(tabClusti,file="tabClusters.csv",append=TRUE,row.names=FALSE,sep=";")  
    }
    
    
    LE_ID_clust=data.frame(LE_ID=LE_ID,clust=metierByLogeventL7)
    print(" --- end of step 3 ---")
    print(Sys.time()-t1)

    return(list(LE_ID_clust=LE_ID_clust, clusters=clusters, sizeClusters=sizeClusters,
     betweenVarClassifOnTot=betweenVarClassifOnTot, mProfilSample=mProfilSample, 
     nbClust=nbClust, summaryClusters=summaryClusters, testValues=resval, 
     testValuesSpecies=target$tabnomespcib, percLogevents=percLogevents,
     descClusters=clusterDesc2, tabClusters=tabClusters,
     targetSpecies=targetSpeciesByCluster))

  }   else



                          

########################################################################################################################################   KMEANS

  if(methMetier=="kmeans"){
    # Calculation of optimal k thanks to within variance
    nbLog=nrow(datLog)

    varintra=numeric()
    for (k in 2:15){
      clustersKmeans=kmeans(datLog, k, iter.max=60, nstart=10)
      varintra[k]=1/nbLog*sum(clustersKmeans$withinss)
    }
    nbClust=which(scree(varintra)[,"epsilon"]<0)[3]
    

    png(paste(analysisName,"Within variance of the classification.png",sep="_"), width = 1200, height = 800)
    plot(varintra,main="Within clusters variance",xlab="Number of clusters",ylab="Within Variance")
    dev.off()

    Store(objects())
    gc(reset=TRUE)

    # KMEANS with k optimal
    clusters=kmeans(datLog, nbClust, iter.max=60, nstart=10, algorithm="Hartigan-Wong")


    # Within and between variance of clusters and classification
    centerOfGravityClassif=numeric()
    withinVarClusters=numeric()
    sizeClusti=numeric()
    centerOfGravityDatLog=colMeans(datLog)
    centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)
    for(k in 1:nbClust){  # Within variance by cluster

      clusti=datLog[which(clusters$cluster==k),]
      if(length(which(clusters$cluster==k))==1)  centerOfGravityClusti=clusti
      else centerOfGravityClusti=colMeans(clusti)
      centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
      sizeClusti[k]=length(which(clusters$cluster==k))
      if(length(which(clusters$cluster==k))==1)  withinVarClusters[k]=0
      else withinVarClusters[k]=sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))

    }
    # Between variance
    classifBetweenVar=(1/nbLog)*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
    # Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
    withinVarClusterOnTot=(1/nbLog)*sum(withinVarClusters)/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
    betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100


    # Compute the test-values for species
    resval=test.values(clusters$cluster,datSpecies)
    # Determine the target species
    target=targetspecies(resval)
    rownames(target$tabnomespcib)=paste("Cluster",1:nbClust)
    
    
    # Compute the percentage of logevents catching each species by cluster
    mainSpecies=colnames(datSpecies)
    percLogevents=matrix(0,ncol=length(mainSpecies),nrow=nbClust,dimnames=list(paste("Cluster ",1:nbClust,sep=""),mainSpecies))
    for(i in 1:nbClust){
      percLogevents[i,]=round(sapply(mainSpecies,function(x) (clusters$size[i]-length(which(datSpecies[clusters$cluster==i,x]==0)))/clusters$size[i]*100),digits=1)
    }
    

    # Projections on the first factorial plans
    png(paste(analysisName,"Projections.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
    if(dim(datLog)[2]>2) {
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
    plot(datLog[,1], datLog[,3], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
    if(dim(datLog)[2]>3) {
    plot(datLog[,1], datLog[,4], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
    plot(datLog[,2], datLog[,4], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
    plot(datLog[,3], datLog[,4], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="Projection of Kmeans classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
    }}
    par(op)
    dev.off()

    # For a paper
    X11(5,5)
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="", xlab="Axis 1", ylab="Axis 2")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_1_2_Kmeans',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="", xlab="Axis 1", ylab="Axis 2")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_1_2_Kmeans',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(clusters$size))[as.numeric(clusters$cluster)], main="", xlab="Axis 2", ylab="Axis 3")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_2_3_Kmeans',sep="_"), type='png')
    dev.off()


    # Catch profile of the dataset
    meanprofile=colMeans(datSpecies)
    png(paste(analysisName,"Catch profile of the dataset.png",sep="_"), width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Catch profile of the dataset", xlab="Species", ylab="Percentage of catch")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    dev.off()


    # Catch profile by cluster
    nbSpec=ncol(datSpecies)
    summaryClusters=array(0,dim=c(6,nbSpec,nbClust))
    dimnames(summaryClusters)[[1]]=c("Min.","1st Qu.","Median", "Mean", "3rd Qu.", "Max.")
    dimnames(summaryClusters)[[2]]=names(meanprofile)
    dimnames(summaryClusters)[[3]]=paste("Cluster",1:nbClust)
    for(i in 1:nbClust){
      if(clusters$size[i]==1){
        summaryClusters[,,i]=apply(t(as.matrix(datSpecies[which(clusters$cluster==i),])),2,
          function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
      }else{
        summaryClusters[,,i]=apply(datSpecies[which(clusters$cluster==i),],2,
          function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
      }
    }
    # Species names for catch profile plots
    nameSpPlot=character()
    catchMeanThreshold=2
    for(i in 1:nbClust){
      namSpi=names(meanprofile[which(t(summaryClusters["Mean",,i])>catchMeanThreshold)])
      numSpi=which(t(summaryClusters["Mean",,i])>catchMeanThreshold)
      nameSpPloti=rep("",nbSpec)
      nameSpPloti[numSpi]=namSpi
      nameSpPlot=rbind(nameSpPlot,nameSpPloti)
    }
    png(paste(analysisName,"Catch profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(t(summaryClusters["Mean",,i]), names.arg=nameSpPlot[i,], xlab="Species", ylab="Percentage of catch", col="gray")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Catch profile by cluster","\n","\n",sep=""))
    dev.off()


    # For a paper : levelplot
    X11(4,4)
    mat <- t(summaryClusters["Mean",,])
    #rownames(mat) <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
    #                          "XVI","XVII","XVIII","XIX","XX")[1:nrow(mat)]
    rownames(mat) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")[1:nrow(mat)]
    sp <- apply(mat,2, function(x) length(which(x>20)))
    sp2=sp[which(sp>=1)]
    #colnames(mat)[sp==0] <- ""
    mat2=mat[,which(colnames(mat)%in%names(sp2))]
    cc <- colorRampPalette(c("white", "black"),space = "rgb", interpolate="spline")
    print(levelplot(mat2, cuts=4, aspect=1, xlab="", ylab="", col.regions=cc(5), at=c(0,20,40,60,80,100), scales=list(cex=0.7), colorkey=list(space="right", at=c(0,20,40,60,80,100), width=1.1)))
    savePlot(filename=paste(analysisName,'mean_profile_by_cluster_levelplot',sep="_"), type='png')
    dev.off()


    # Standard deviation profile by cluster
    sdprofil=matrix(0,nrow=nbClust,ncol=nbSpec)
    namSdPlot=character()
    SdThreshold=2
    for(i in 1:nbClust){
      if(length(which(clusters$cluster==i))==1){ sdprofilclusti=rep(0,nbSpec)
      }else{sdprofilclusti=apply(datSpecies[which(clusters$cluster==i),],2,sd)}
      namSDi=names(which(sdprofilclusti>SdThreshold))
      numSDi=which(sdprofilclusti>SdThreshold)
      namSdPloti=rep("",nbSpec)
      namSdPloti[numSDi]=namSDi
      sdprofil[i,]=sdprofilclusti
      namSdPlot=rbind(namSdPlot,namSdPloti)
    }
    rownames(sdprofil) <- 1:nrow(sdprofil)
    png(paste(analysisName,"Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(sdprofil[i,], names.arg=namSdPlot[i,], xlab="Species", ylab="Percentage of catch")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Standard deviation profile by cluster","\n","\n",sep=""))
    dev.off()


    # Number of Logevents by cluster
    x=c(1:nbClust)
    png(paste(analysisName,"Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
    coord=barplot(clusters$size, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    barplot(clusters$size, names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
    text(coord,clusters$size+400,clusters$size,font=2,xpd=NA)
    dev.off()


    # Profile of test-values by cluster
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

    png(paste(analysisName,"Profile of test-values by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,], cex.names=1, names.arg=nameTargetPlot[i,], xlab="Species", ylab="Test-value")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Profile of test-values by cluster","\n","\n",sep=""))
    dev.off()




    # Descriptive tables of the clusters
    clusterDesc=matrix(0,nrow=9,ncol=nbClust)
    for(i in 1:nbClust){
      clusterDesc[,i]=c(i,
                        length(which(cumsum(t(summaryClusters["Mean",,i])[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)])<50))+1,
                        length(which(cumsum(t(summaryClusters["Mean",,i])[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)])<90))+1,
                        length(which(t(summaryClusters["Median",,i])>50)),
                        length(which(resval[,i]>1.96)),
                        length(which(resval[,i]>3.29)),
                        length(which(apply(datSpecies,2,function (x) (clusters$size[i]-length(which(x[clusters$cluster==i]==0)))/clusters$size[i]*100)>50)),
                        length(which(apply(datSpecies,2,function (x) (clusters$size[i]-length(which(x[clusters$cluster==i]==0)))/clusters$size[i]*100)>90)),
                        clusters$size[i])
    }
    rownames(clusterDesc)=c("Number of species",
                            "to have 50% of catch", "to have 90% of catch",
                            "with a median higher than 50",
                            "with a test-value > 1.96", "with a test-value > 3.29",
                            "catch in 50% of the logevents", "catch in 90% of the logevents",
                            "Clusters size")
    colnames(clusterDesc)=1:nbClust
    clusterDesc2=as.data.frame(clusterDesc)


    # Summary tables of the clusters
    namesSpecies=matrix(NA,nrow=nbClust,ncol=10)
    namesCapt=matrix(NA,nrow=nbClust,ncol=5)
    nbSpeciesCatch = min(5,dim(t(summaryClusters["Mean",,]))[[2]])
    namesTarget=matrix(NA,nrow=nbClust,ncol=5)
    nbSpeciesVT = min(5,dim(target$tabnomespcib)[[2]])
    tabLibname=matrix(NA,nrow=nbClust,ncol=10)
    listLibname=list()

    for(i in 1:nbClust){
      namesCapt[i,]=colnames(t(summaryClusters["Mean",,i]))[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)][1:nbSpeciesCatch]
      a=as.data.frame(t(summaryClusters["Mean",target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])],i]))
      colnames(a)= target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])]
      if(length(a)!=0){
        namesTarget[i,1:length(target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])])]=colnames(a[order(a,decreasing=TRUE)])
      }
      namesSpecies[i,1:length(union(namesCapt[i,],namesTarget[i,]))]=union(namesCapt[i,],namesTarget[i,])
    }
    
    for(i in 1:nbClust){
        listLibname[[i]]=lapply(as.list(namesSpecies[i,]), function(x) if(length(which(correspLevel7to5[,"X3A_CODE"]==x))==0) "NA" 
                                                    else correspLevel7to5[which(correspLevel7to5[,"X3A_CODE"]==x),"French_name"])
        tabLibname[i,]=unlist(lapply(listLibname[[i]], function(x) as.character(unlist(x))))
    }
        
    tabPropCatch=matrix(NA,nrow=nbClust,ncol=10)
    tabTestVal=matrix(NA,nrow=nbClust,ncol=10)
    tabPropLog=matrix(NA,nrow=nbClust,ncol=10)

    for(i in 1:nbClust){
      print("-----------------------------------------------------------------")
      print(paste("Cluster",i))
      propCatch=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)>=0.1)]
      tabPropCatch[i,1:length(propCatch)]=propCatch
      print(propCatch)
      testVal=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) resval[x,i]),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)>=0.1)]
      tabTestVal[i,1:length(testVal)]=testVal
      print(testVal)
      propLog=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) (clusters$size[i]-length(which(datSpecies[clusters$cluster==i,x]==0)))/clusters$size[i]*100),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) (clusters$size[i]-length(which(datSpecies[clusters$cluster==i,x]==0)))/clusters$size[i]*100),digits=1)>=0.1)]
      tabPropLog[i,1:length(propLog)]=propLog
      print(propLog)
    }

    tabClusters=array(0,dim=c(10,5,nbClust))
    dimnames(tabClusters)[[2]]=c("Libname","FAO","Test-value","% Catch","% Logevents")
    dimnames(tabClusters)[[3]]=paste("Cluster",1:nbClust)
    for(i in 1:nbClust){
      tabClusters[,,i]=cbind(tabLibname[i,],namesSpecies[i,],tabTestVal[i,],tabPropCatch[i,],tabPropLog[i,])
    }

    sizeTabClusters=numeric()
    for(i in 1:nbClust){
      sizeTabClusters[i]=min(length(namesSpecies[i,!is.na(namesSpecies[i,])]),length(tabPropCatch[i,!is.na(tabPropCatch[i,])]),length(tabTestVal[i,!is.na(tabTestVal[i,])]),length(tabPropLog[i,!is.na(tabPropLog[i,])]))
    }


    # Target Species
    # Intersection of species from tabClusters having : - % Cumulated Catch > thresholdCatch
    #                                                   - Test-value > thresholdTestValue
    #                                                   - % Logevents > thresholdLogevents 
    thresholdCatch=75
    thresholdTestValue=3
    thresholdLogevents=30
    
    sppCumCatch=list()
    sppTestValue=list()
    sppLogevents=list()
    targetSpeciesByCluster=list()
    
    for (i in 1:nbClust){
      percCatchCum=cumsum(as.numeric(tabClusters[,"% Catch",i]))
      nbSpSel=length(which(percCatchCum<thresholdCatch))+1
      sppCumCatch[[i]]=tabClusters[1:nbSpSel,"FAO",i]
      
      sppTestValue[[i]]=tabClusters[which(as.numeric(tabClusters[,"Test-value",i])>thresholdTestValue),"FAO",i]
      
      sppLogevents[[i]]=tabClusters[which(as.numeric(tabClusters[,"% Logevents",i])>thresholdLogevents),"FAO",i]
    
      targetSpeciesByCluster[[i]]=intersect(sppCumCatch[[i]],sppTestValue[[i]])
      targetSpeciesByCluster[[i]]=intersect(targetSpeciesByCluster[[i]],sppLogevents[[i]])
    }

    # List of metiers (level 7)
    listMetiersL7=list()
    for (i in 1:nbClust){
      metiersClusteri=targetSpeciesByCluster[[i]]
      metiersClusteri=as.character(unique(unlist(metiersClusteri)))
      metiersClusteri=paste(unlist(strsplit(metiersClusteri," ")),collapse=" ")
      listMetiersL7[[i]]=metiersClusteri
    }
  
    # Metier (level 7) of each logevent
    metierByLogeventL7=unlist(sapply(clusters$cluster,function(x) listMetiersL7[[x]]))



    # Create csv tables
    write.table(clusterDesc2,file="descClusters.csv",col.names=FALSE,sep=";")
    
    dfClust=data.frame()
    dfClust=paste("Clust ",1:nbClust,sep="")
    for(i in 1:nbClust){
      write.table(dfClust[i],file="tabClusters.csv",append=TRUE,row.names=FALSE,sep=";")
      tabClusti=as.data.frame(tabClusters[1:sizeTabClusters[i],,i])
      write.table(tabClusti,file="tabClusters.csv",append=TRUE,row.names=FALSE,sep=";")  
    }
    
    
    LE_ID_clust=data.frame(LE_ID=LE_ID,clust=metierByLogeventL7)
    print(" --- end of step 3 ---")
    print(Sys.time()-t1)

    return(list(LE_ID_clust=LE_ID_clust, clusters=clusters, 
    betweenVarClassifOnTot=betweenVarClassifOnTot, nbClust=nbClust, 
    summaryClusters=summaryClusters, testValues=resval, 
    testValuesSpecies=target$tabnomespcib, percLogevents=percLogevents,
    descClusters=clusterDesc2, tabClusters=tabClusters,
    targetSpecies=targetSpeciesByCluster))

  } else


                      


########################################################################################################################################   PAM

  if(methMetier=="pam"){

    # Calculation of optimal k thanks to the silhouette (second maximum)
    nbLog=nrow(datLog)

    clustersPam.silcoeff=numeric()
    clustersPam.silcoeff[1]=0
    clustersPam.silcoeff[2]=0
    clustersPam.silcoeff[3]=0
    
    k=2
    compMax=1
    repeat{
      k=k+2
      print(k)
      clustersPam=pam(datLog,k)
      clustersPam.silcoeff[k]=clustersPam$silinfo$avg.width
      clustersPam=pam(datLog,k+1)
      clustersPam.silcoeff[k+1]=clustersPam$silinfo$avg.width
      if((clustersPam.silcoeff[k-2]<clustersPam.silcoeff[k-1] & clustersPam.silcoeff[k-1]>clustersPam.silcoeff[k]) & compMax<=2){
        if(compMax==2){
          nbClust=k-1
          print(paste("2e max =",k-1))
          print(paste("nbClust =",nbClust))
          break
        } else {
          compMax=compMax+1
          print(paste("compMax1 =",compMax))
          print(paste("1er max =",k-1))
        }
      }
      if((clustersPam.silcoeff[k-1]<clustersPam.silcoeff[k] & clustersPam.silcoeff[k]>clustersPam.silcoeff[k+1]) & compMax<=2){
        if(compMax==2){
          nbClust=k
          print(paste("2e max =",k))
          print(paste("nbClust =",nbClust))
          break
        } else {
          compMax=compMax+1
          print(paste("compMax2 =",compMax))
          print(paste("1er max =",k))
        }
      }
    Store(objects())
    gc(reset=TRUE)
    }

    png(paste(analysisName,"Silhouette of the classification.png",sep="_"), width = 1200, height = 800)
    plot(clustersPam.silcoeff, main="Silhouette of the classification", xlab="Number of clusters", ylab="Silhouette") # k optimal corresponds to second maximum of silhouette's coefficients
    dev.off()

    Store(objects())
    gc(reset=TRUE)

    cat("PamSilCoeff",clustersPam.silcoeff,"\n")

    # PAM with optimal k
    clusters=pam(datLog,nbClust)
    summary(clusters)


    # Within and between variance of clusters and classification
    centerOfGravityClassif=numeric()
    withinVarClusters=numeric()
    sizeClusti=numeric()
    centerOfGravityDatLog=colMeans(datLog)
    centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)
    for(k in 1:nbClust){  # Within variance by cluster

      clusti=datLog[which(clusters$clustering==k),]
      if(length(which(clusters$clustering==k))==1)  centerOfGravityClusti=clusti
      else centerOfGravityClusti=colMeans(clusti)
      centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
      sizeClusti[k]=length(which(clusters$clustering==k))
      if(length(which(clusters$clustering==k))==1)  withinVarClusters[k]=0
      else withinVarClusters[k]=sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))

    }
    # Between variance
    classifBetweenVar=(1/nbLog)*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
    # Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
    withinVarClusterOnTot=(1/nbLog)*sum(withinVarClusters)/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
    betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100


    # Compute the test-values for species
    resval=test.values(clusters$clustering,datSpecies)
    # Determine the target species
    target=targetspecies(resval)
    nbClust=length(clusters$id.med)
    rownames(target$tabnomespcib)=paste("Cluster",1:nbClust)
    
    
    # Compute the percentage of logevents catching each species by cluster
    mainSpecies=colnames(datSpecies)
    percLogevents=matrix(0,ncol=length(mainSpecies),nrow=nbClust,dimnames=list(paste("Cluster ",1:nbClust,sep=""),mainSpecies))
    for(i in 1:nbClust){
      percLogevents[i,]=round(sapply(mainSpecies,function(x) (clusters$clusinfo[i,1]-length(which(datSpecies[clusters$clustering==i,x]==0)))/clusters$clusinfo[i,1]*100),digits=1)
    }
    

    # Projections on the first factorial plans
    png(paste(analysisName,"Projections.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
    if(dim(datLog)[2]>2) {
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
    plot(datLog[,1], datLog[,3], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
    if(dim(datLog)[2]>3) {
    plot(datLog[,1], datLog[,4], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
    plot(datLog[,2], datLog[,4], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
    plot(datLog[,3], datLog[,4], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="Projection of PAM classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
    }}
    par(op)
    dev.off()


    # For a paper
    X11(5,5)
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="", xlab="Axis 1", ylab="Axis 2")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_1_2_PAM',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="", xlab="Axis 1", ylab="Axis 2")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_1_2_PAM',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(clusters$id.med))[as.numeric(clusters$clustering)], main="", xlab="Axis 2", ylab="Axis 3")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_2_3_PAM',sep="_"), type='png')
    dev.off()

    
    # Catch profile of the dataset
    meanprofile=colMeans(datSpecies)
    png(paste(analysisName,"Catch profile of the dataset.png",sep="_"), width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Catch profile of the dataset", xlab="Species", ylab="Percentage of catch")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    dev.off()


    # Catch profile by cluster
    nbSpec=ncol(datSpecies)
    summaryClusters=array(0,dim=c(6,nbSpec,nbClust))
    dimnames(summaryClusters)[[1]]=c("Min.","1st Qu.","Median", "Mean", "3rd Qu.", "Max.")
    dimnames(summaryClusters)[[2]]=names(meanprofile)
    dimnames(summaryClusters)[[3]]=paste("Cluster",1:nbClust)
    for(i in 1:nbClust){
      if(clusters$clusinfo[i,1]==1){
        summaryClusters[,,i]=apply(t(as.matrix(datSpecies[which(clusters$clustering==i),])),2,
          function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
      }else{
        summaryClusters[,,i]=apply(datSpecies[which(clusters$clustering==i),],2,
          function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
      }
    }
    # Species names for catch profile plots
    nameSpPlot=character()
    catchMeanThreshold=2
    for(i in 1:nbClust){
      namSpi=names(meanprofile[which(t(summaryClusters["Mean",,i])>catchMeanThreshold)])
      numSpi=which(t(summaryClusters["Mean",,i])>catchMeanThreshold)
      nameSpPloti=rep("",nbSpec)
      nameSpPloti[numSpi]=namSpi
      nameSpPlot=rbind(nameSpPlot,nameSpPloti)
    }
    png(paste(analysisName,"Catch profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(t(summaryClusters["Mean",,i]), names.arg=nameSpPlot[i,], xlab="Species", ylab="Percentage of catch", col="gray")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Catch profile by cluster","\n","\n",sep=""))
    dev.off()


    # For a paper : levelplot
    X11(4,4)
    mat <- t(summaryClusters["Mean",,])
    #rownames(mat) <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
    #                          "XVI","XVII","XVIII","XIX","XX")[1:nrow(mat)]
    rownames(mat) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")[1:nrow(mat)]
    sp <- apply(mat,2, function(x) length(which(x>20)))
    sp2=sp[which(sp>=1)]
    #colnames(mat)[sp==0] <- ""
    mat2=mat[,which(colnames(mat)%in%names(sp2))]
    cc <- colorRampPalette(c("white", "black"),space = "rgb", interpolate="spline")
    print(levelplot(mat2, cuts=4, aspect=1, xlab="", ylab="", col.regions=cc(5), at=c(0,20,40,60,80,100), scales=list(cex=0.7), colorkey=list(space="right", at=c(0,20,40,60,80,100), width=1.1)))
    savePlot(filename=paste(analysisName,'mean_profile_by_cluster_levelplot',sep="_"), type='png')
    dev.off()


    # Standard deviation profile by cluster
    sdprofil=matrix(0,nrow=nbClust,ncol=nbSpec)
    namSdPlot=character()
    SdThreshold=2
    for(i in 1:nbClust){
      if(length(which(clusters$clustering==i))==1){ sdprofilclusti=rep(0,nbSpec)
      }else{sdprofilclusti=apply(datSpecies[which(clusters$clustering==i),],2,sd)}
      namSDi=names(which(sdprofilclusti>SdThreshold))
      numSDi=which(sdprofilclusti>SdThreshold)
      namSdPloti=rep("",nbSpec)
      namSdPloti[numSDi]=namSDi
      sdprofil[i,]=sdprofilclusti
      namSdPlot=rbind(namSdPlot,namSdPloti)
    }
    rownames(sdprofil) <- 1:nrow(sdprofil)
    png(paste(analysisName,"Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(sdprofil[i,], names.arg=namSdPlot[i,], xlab="Species", ylab="Percentage of catch")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Standard deviation profile by cluster","\n","\n",sep=""))
    dev.off()


    # Number of Logevents by cluster
    x=c(1:nbClust)
    png(paste(analysisName,"Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
    coord=barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
    text(coord,clusters$clusinfo[,1]+200,clusters$clusinfo[,1],font=2,xpd=NA)
    dev.off()


    # Profile of test-values by cluster
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

    png(paste(analysisName,"Profile of test-values by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,],names.arg=nameTargetPlot[i,], xlab="Species", ylab="Test-value")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Profile of test-values by cluster","\n","\n",sep=""))
    dev.off()



    # Descriptive tables of the clusters
    clusterDesc=matrix(0,nrow=9,ncol=nbClust)
    for(i in 1:nbClust){
      clusterDesc[,i]=c(i,
                        length(which(cumsum(t(summaryClusters["Mean",,i])[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)])<50))+1,
                        length(which(cumsum(t(summaryClusters["Mean",,i])[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)])<90))+1,
                        length(which(t(summaryClusters["Median",,i])>50)),
                        length(which(resval[,i]>1.96)),
                        length(which(resval[,i]>3.29)),
                        length(which(apply(datSpecies,2,function (x) (clusters$clusinfo[i,1]-length(which(x[clusters$clustering==i]==0)))/clusters$clusinfo[i,1]*100)>50)),
                        length(which(apply(datSpecies,2,function (x) (clusters$clusinfo[i,1]-length(which(x[clusters$clustering==i]==0)))/clusters$clusinfo[i,1]*100)>90)),
                        clusters$clusinfo[i,1])
    }
    rownames(clusterDesc)=c("Number of species",
                            "to have 50% of catch", "to have 90% of catch",
                            "with a median higher than 50",
                            "with a test-value > 1.96", "with a test-value > 3.29",
                            "catch in 50% of the logevents", "catch in 90% of the logevents",
                            "Clusters size")
    colnames(clusterDesc)=1:nbClust
    clusterDesc2=as.data.frame(clusterDesc)


    # Summary tables of the clusters
    namesSpecies=matrix(NA,nrow=nbClust,ncol=10)
    namesCapt=matrix(NA,nrow=nbClust,ncol=5)
    nbSpeciesCatch = min(5,dim(t(summaryClusters["Mean",,]))[[2]])
    namesTarget=matrix(NA,nrow=nbClust,ncol=5)
    nbSpeciesVT = min(5,dim(target$tabnomespcib)[[2]])
    tabLibname=matrix(NA,nrow=nbClust,ncol=10)
    listLibname=list()

    for(i in 1:nbClust){
      namesCapt[i,]=colnames(t(summaryClusters["Mean",,i]))[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)][1:nbSpeciesCatch]
      a=as.data.frame(t(summaryClusters["Mean",target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])],i]))
      colnames(a)= target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])]
      if(length(a)!=0){
        namesTarget[i,1:length(target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])])]=colnames(a[order(a,decreasing=TRUE)])
      }
      namesSpecies[i,1:length(union(namesCapt[i,],namesTarget[i,]))]=union(namesCapt[i,],namesTarget[i,])
    }
    
    for(i in 1:nbClust){
        listLibname[[i]]=lapply(as.list(namesSpecies[i,]), function(x) if(length(which(correspLevel7to5[,"X3A_CODE"]==x))==0) "NA" 
                                                    else correspLevel7to5[which(correspLevel7to5[,"X3A_CODE"]==x),"French_name"])
        tabLibname[i,]=unlist(lapply(listLibname[[i]], function(x) as.character(unlist(x))))
    }

    tabPropCatch=matrix(NA,nrow=nbClust,ncol=10)
    tabTestVal=matrix(NA,nrow=nbClust,ncol=10)
    tabPropLog=matrix(NA,nrow=nbClust,ncol=10)

    for(i in 1:nbClust){
      print("-----------------------------------------------------------------")
      print(paste("Cluster",i))
      propCatch=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)>=0.1)]
      tabPropCatch[i,1:length(propCatch)]=propCatch
      print(propCatch)
      testVal=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) resval[x,i]),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)>=0.1)]
      tabTestVal[i,1:length(testVal)]=testVal
      print(testVal)
      propLog=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) (clusters$clusinfo[i,1]-length(which(datSpecies[clusters$clustering==i,x]==0)))/clusters$clusinfo[i,1]*100),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) (clusters$clusinfo[i,1]-length(which(datSpecies[clusters$clustering==i,x]==0)))/clusters$clusinfo[i,1]*100),digits=1)>=0.1)]
      tabPropLog[i,1:length(propLog)]=propLog
      print(propLog)
    }

    tabClusters=array(0,dim=c(10,5,nbClust))
    dimnames(tabClusters)[[2]]=c("Libname","FAO","Test-value","% Catch","% Logevents")
    dimnames(tabClusters)[[3]]=paste("Cluster",1:nbClust)
    for(i in 1:nbClust){
      tabClusters[,,i]=cbind(tabLibname[i,],namesSpecies[i,],tabTestVal[i,],tabPropCatch[i,],tabPropLog[i,])
    }

    sizeTabClusters=numeric()
    for(i in 1:nbClust){
      sizeTabClusters[i]=min(length(namesSpecies[i,!is.na(namesSpecies[i,])]),length(tabPropCatch[i,!is.na(tabPropCatch[i,])]),length(tabTestVal[i,!is.na(tabTestVal[i,])]),length(tabPropLog[i,!is.na(tabPropLog[i,])]))
    }


    # Target Species
    # Intersection of species from tabClusters having : - % Cumulated Catch > thresholdCatch
    #                                                   - Test-value > thresholdTestValue
    #                                                   - % Logevents > thresholdLogevents 
    thresholdCatch=75
    thresholdTestValue=3
    thresholdLogevents=30
    
    sppCumCatch=list()
    sppTestValue=list()
    sppLogevents=list()
    targetSpeciesByCluster=list()
    
    for (i in 1:nbClust){
      percCatchCum=cumsum(as.numeric(tabClusters[,"% Catch",i]))
      nbSpSel=length(which(percCatchCum<thresholdCatch))+1
      sppCumCatch[[i]]=tabClusters[1:nbSpSel,"FAO",i]
      
      sppTestValue[[i]]=tabClusters[which(as.numeric(tabClusters[,"Test-value",i])>thresholdTestValue),"FAO",i]
      
      sppLogevents[[i]]=tabClusters[which(as.numeric(tabClusters[,"% Logevents",i])>thresholdLogevents),"FAO",i]
    
      targetSpeciesByCluster[[i]]=intersect(sppCumCatch[[i]],sppTestValue[[i]])
      targetSpeciesByCluster[[i]]=intersect(targetSpeciesByCluster[[i]],sppLogevents[[i]])
    }

    # List of metiers (level 7)
    listMetiersL7=list()
    for (i in 1:nbClust){
      metiersClusteri=targetSpeciesByCluster[[i]]
      metiersClusteri=as.character(unique(unlist(metiersClusteri)))
      metiersClusteri=paste(unlist(strsplit(metiersClusteri," ")),collapse=" ")
      listMetiersL7[[i]]=metiersClusteri
    }
  
    # Metier (level 7) of each logevent
    metierByLogeventL7=unlist(sapply(clusters$clustering,function(x) listMetiersL7[[x]]))


    # Create csv tables
    write.table(clusterDesc2,file="descClusters.csv",col.names=FALSE,sep=";")
    
    dfClust=data.frame()
    dfClust=paste("Clust ",1:nbClust,sep="")
    for(i in 1:nbClust){
      write.table(dfClust[i],file="tabClusters.csv",append=TRUE,row.names=FALSE,sep=";")
      tabClusti=as.data.frame(tabClusters[1:sizeTabClusters[i],,i])
      write.table(tabClusti,file="tabClusters.csv",append=TRUE,row.names=FALSE,sep=";")  
    }


    LE_ID_clust=data.frame(LE_ID=LE_ID,clust=metierByLogeventL7)
    print(" --- end of step 3 ---")
    print(Sys.time()-t1)

    return(list(LE_ID_clust=LE_ID_clust, clusters=clusters, 
    betweenVarClassifOnTot=betweenVarClassifOnTot, nbClust=nbClust, 
    summaryClusters=summaryClusters, testValues=resval, 
    testValuesSpecies=target$tabnomespcib, percLogevents=percLogevents,
    descClusters=clusterDesc2, tabClusters=tabClusters,
    targetSpecies=targetSpeciesByCluster))

  } else


         


########################################################################################################################################   CLARA

  if(methMetier=="clara"){
    nbLog=nrow(datLog)
    propSample=0.1

    # Calculation of optimal k thanks to the silhouette (second maximum)
    clustersClara.silcoeff=numeric()
    clustersClara.silcoeff[1]=0
    clustersClara.silcoeff[2]=0
    clustersClara.silcoeff[3]=0
    k=2
    compMax=1
    repeat{
      k=k+2
      print(k)
      clustersClara=clara(datLog, k, metric=param1, stand=FALSE, samples=5, sampsize=min(nbLog,round(propSample*nbLog+10*k)))
      clustersClara.silcoeff[k]=clustersClara$silinfo$avg.width
      clustersClara=clara(datLog, k+1, metric=param1, stand=FALSE, samples=5, sampsize=min(nbLog,round(propSample*nbLog+10*(k+1))))
      clustersClara.silcoeff[k+1]=clustersClara$silinfo$avg.width
      if((clustersClara.silcoeff[k-2]<clustersClara.silcoeff[k-1] & clustersClara.silcoeff[k-1]>clustersClara.silcoeff[k]) & compMax<=2){
        if(compMax==2){
          nbClust=k-1
          print(paste("2e max =",k-1))
          print(paste("nbClust =",nbClust))
          break
        } else {
          compMax=compMax+1
          print(paste("compMax1 =",compMax))
          print(paste("1er max =",k-1))
        }
      }
      if((clustersClara.silcoeff[k-1]<clustersClara.silcoeff[k] & clustersClara.silcoeff[k]>clustersClara.silcoeff[k+1]) & compMax<=2){
        if(compMax==2){
          nbClust=k
          print(paste("2e max =",k))
          print(paste("nbClust =",nbClust))
          break
        } else {
          compMax=compMax+1
          print(paste("compMax2 =",compMax))
          print(paste("1er max =",k))
        }
      }
    Store(objects())
    gc(reset=TRUE)
    }
    

    png(paste(analysisName,"Silhouette of the classification.png",sep="_"), width = 1200, height = 800)
    plot(clustersClara.silcoeff, main="Silhouette of the classification", xlab="Number of clusters", ylab="Silhouette")               # k optimal corresponds to maximum of silhouette's coefficients
    dev.off()

    Store(objects())
    gc(reset=TRUE)

    cat("ClaraSilCoeff",clustersClara.silcoeff,"\n")


    # CLARA with optimal k
    clusters=clara(datLog, nbClust, metric=param1, stand=FALSE, samples=5, sampsize=min(nbLog,round(propSample*nbLog+10*nbClust)))  # CLARA with optimal k
    summary(clusters)


    # Within and between variance of clusters and classification
    centerOfGravityClassif=numeric()
    withinVarClusters=numeric()
    sizeClusti=numeric()
    centerOfGravityDatLog=colMeans(datLog)
    centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityDatLog)
    for(k in 1:nbClust){  # Within variance by cluster

      clusti=datLog[which(clusters$clustering==k),]
      if(length(which(clusters$clustering==k))==1)  centerOfGravityClusti=clusti
      else centerOfGravityClusti=colMeans(clusti)
      centerOfGravityClassif=rbind(centerOfGravityClassif,centerOfGravityClusti)
      sizeClusti[k]=length(which(clusters$clustering==k))
      if(length(which(clusters$clustering==k))==1)  withinVarClusters[k]=0
      else withinVarClusters[k]=sum(apply(clusti,1,function(x) withinVar(x,centerOfGravityClusti)))
      
    }
    # Between variance
    classifBetweenVar=(1/nbLog)*sum(sizeClusti*((dist(centerOfGravityClassif)[1:nbClust])^2))
    # Within variance of clusters on totale variance (pourcent) and between variance on totale variance of classification
    withinVarClusterOnTot=(1/nbLog)*sum(withinVarClusters)/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100
    betweenVarClassifOnTot=classifBetweenVar/(classifBetweenVar+(1/nbLog)*sum(withinVarClusters))*100


    # Compute the test-values for species
    resval=test.values(clusters$cluster,datSpecies)
    # Determine the target species
    target=targetspecies(resval)
    rownames(target$tabnomespcib)=paste("Cluster",1:nbClust)
    
    
    # Compute the percentage of logevents catching each species by cluster
    mainSpecies=colnames(datSpecies)
    percLogevents=matrix(0,ncol=length(mainSpecies),nrow=nbClust,dimnames=list(paste("Cluster ",1:nbClust,sep=""),mainSpecies))
    for(i in 1:nbClust){
      percLogevents[i,]=round(sapply(mainSpecies,function(x) (clusters$clusinfo[i,1]-length(which(datSpecies[clusters$clustering==i,x]==0)))/clusters$clusinfo[i,1]*100),digits=1)
    }
    
    
    # Projections on the first factorial plans
    png(paste(analysisName,"Projections.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 1-2", xlab="axis 1", ylab="axis 2")
    if(dim(datLog)[2]>2) {
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 2-3", xlab="axis 2", ylab="axis 3")
    plot(datLog[,1], datLog[,3], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 1-3", xlab="axis 1", ylab="axis 3")
    if(dim(datLog)[2]>3) {
    plot(datLog[,1], datLog[,4], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 1-4", xlab="axis 1", ylab="axis 4")
    plot(datLog[,2], datLog[,4], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 2-4", xlab="axis 2", ylab="axis 4")
    plot(datLog[,3], datLog[,4], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="Projection of CLARA classification on the factorial plan 3-4", xlab="axis 3", ylab="axis 4")
    }}
    par(op)
    dev.off()


    # For a paper
    X11(5,5)
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="", xlab="Axis 1", ylab="Axis 2")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_1_2_CLARA',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot(datLog[,1], datLog[,2], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="", xlab="Axis 1", ylab="Axis 2")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_1_2_CLARA',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot(datLog[,2], datLog[,3], pch=21, bg=rainbow(length(clusters$i.med))[as.numeric(clusters$clustering)], main="", xlab="Axis 2", ylab="Axis 3")
    abline(h=0, lty=2) ; abline(v=0, lty=2)
    savePlot(filename=paste(analysisName,'projections_2_3_CLARA',sep="_"), type='png')
    dev.off()

    # Catch profile of the dataset
    meanprofile=colMeans(datSpecies)
    png(paste(analysisName,"Catch profile of the dataset.png",sep="_"), width = 1200, height = 800)
    op <- par(las=2)
    barplot(meanprofile, main="Catch profile of the dataset", xlab="Species", ylab="Percentage of catch")
    par(op)
    mtext(paste(nrow(datSpecies)," logevents"), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    dev.off()


    # Catch profile by cluster
    nbSpec=ncol(datSpecies)
    summaryClusters=array(0,dim=c(6,nbSpec,nbClust))
    dimnames(summaryClusters)[[1]]=c("Min.","1st Qu.","Median", "Mean", "3rd Qu.", "Max.")
    dimnames(summaryClusters)[[2]]=names(meanprofile)
    dimnames(summaryClusters)[[3]]=paste("Cluster",1:nbClust)
    for(i in 1:nbClust){
      if(clusters$clusinfo[i,1]==1){
        summaryClusters[,,i]=apply(t(as.matrix(datSpecies[which(clusters$clustering==i),])),2,
          function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
      }else{
        summaryClusters[,,i]=apply(datSpecies[which(clusters$clustering==i),],2,
          function(x) rbind(min(as.vector(x)),quantile(as.vector(x),0.25),quantile(as.vector(x),0.50),mean(as.vector(x)),quantile(as.vector(x),0.75),max(as.vector(x))))
      }
    }
    # Species names for catch profile plots
    nameSpPlot=character()
    catchMeanThreshold=2
    for(i in 1:nbClust){
      namSpi=names(meanprofile[which(t(summaryClusters["Mean",,i])>catchMeanThreshold)])
      numSpi=which(t(summaryClusters["Mean",,i])>catchMeanThreshold)
      nameSpPloti=rep("",nbSpec)
      nameSpPloti[numSpi]=namSpi
      nameSpPlot=rbind(nameSpPlot,nameSpPloti)
    }
    # Plot
    png(paste(analysisName,"Catch profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(t(summaryClusters["Mean",,i]), names.arg=nameSpPlot[i,], xlab="Species", ylab="Percentage of catch", col="gray")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Catch profile by cluster","\n","\n",sep=""))
    dev.off()


    # For a paper : levelplot
    X11(4,4)
    mat <- t(summaryClusters["Mean",,])
    #rownames(mat) <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV",
    #                          "XVI","XVII","XVIII","XIX","XX")[1:nrow(mat)]
    rownames(mat) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")[1:nrow(mat)]
    sp <- apply(mat,2, function(x) length(which(x>20)))
    sp2=sp[which(sp>=1)]
    #colnames(mat)[sp==0] <- ""
    mat2=mat[,which(colnames(mat)%in%names(sp2))]
    cc <- colorRampPalette(c("white", "black"),space = "rgb", interpolate="spline")
    print(levelplot(mat2, cuts=4, aspect=1, xlab="", ylab="", col.regions=cc(5), at=c(0,20,40,60,80,100), scales=list(cex=0.7), colorkey=list(space="right", at=c(0,20,40,60,80,100), width=1.1)))
    savePlot(filename=paste(analysisName,'mean_profile_by_cluster_levelplot',sep="_"), type='png')
    dev.off()

    # OR #
    mat <- t(summaryClusters["Mean",,])
    #mat=mat[,order(colnames(mat),decreasing=TRUE)]     # si on veut mettre les especes par ordre descendant
    rownames(mat) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")[1:nrow(mat)]
    sp <- apply(mat,2, function(x) length(which(x>20)))
    colnames(mat)[sp==0] <- ""
    cc <- colorRampPalette(c("white", "steelblue2", "blue4"),space = "rgb", interpolate="spline")
    png(filename = paste(paste(analysisName,'mean_profile_by_cluster_levelplot_blue',sep="_"),".png",sep=""), width = 400, height = 800)
    print(levelplot(mat, cuts=4, aspect=3, xlab="", ylab="", col.regions=cc(5), at=c(0,20,40,60,80,100), scales=list(cex=0.8), colorkey=list(space="right", at=c(0,20,40,60,80,100))))
    dev.off()
  

    # Standard deviation profile by cluster
    sdprofil=matrix(0,nrow=nbClust,ncol=nbSpec)
    namSdPlot=character()
    SdThreshold=5
    for(i in 1:nbClust){
      if(length(which(clusters$clustering==i))==1){ sdprofilclusti=rep(0,nbSpec)
      }else{sdprofilclusti=apply(datSpecies[which(clusters$clustering==i),],2,sd)}
      namSDi=names(which(sdprofilclusti>SdThreshold))
      numSDi=which(sdprofilclusti>SdThreshold)
      namSdPloti=rep("",nbSpec)
      namSdPloti[numSDi]=namSDi
      sdprofil[i,]=sdprofilclusti
      namSdPlot=rbind(namSdPlot,namSdPloti)
    }
    rownames(sdprofil) <- 1:nrow(sdprofil)
    png(paste(analysisName,"Standard deviation profile by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(sdprofil[i,], names.arg=namSdPlot[i,], xlab="Species", ylab="Percentage of catch")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Standard deviation profile by cluster","\n","\n",sep=""))
    dev.off()


    # Number of Logevents by cluster
    x=c(1:nbClust)
    png(paste(analysisName,"Number of Logevents by cluster.png",sep="_"), width = 1200, height = 800)
    coord=barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents")
    barplot(clusters$clusinfo[,1], names.arg=x, main="Number of Logevents by cluster", xlab="Cluster", ylab="Number of Logevents", col="skyblue")
    text(coord,clusters$clusinfo[,1]+5,clusters$clusinfo[,1],font=2,xpd=NA)
    dev.off()


    # Profile of test-values by cluster
    targetresval=matrix(0,nrow=nbClust,ncol=nbSpec)
    colnames(targetresval)=colnames(datSpecies)
    rownames(targetresval)=1:nbClust
    nameTargetPlot=matrix(NA,nrow=nbClust,ncol=nbSpec)
    for(i in 1:nbClust){
      nomtargeti=as.character(target$tabnomespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      numtargeti=as.numeric(target$tabnumespcib[i,which(!is.na(target$tabnumespcib[i,]))])
      nameTargetPlot[i,numtargeti]=nomtargeti
      targetresval[i,numtargeti]=resval[nomtargeti,i]
    }

    png(paste(analysisName,"Profile of test-values by cluster.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(ceiling(sqrt(nbClust)),round(sqrt(nbClust))))
    for(i in 1:nbClust){
      op2 <- par(las=2)
      barplot(targetresval[i,],names.arg=nameTargetPlot[i,], xlab="Species", ylab="Test-value")
      par(op2)
      mtext(paste("Cluster",i), side=3, outer=FALSE, adj=0.5, line=0.5, col="darkblue")
    }
    par(op)
    title(main=paste("Profile of test-values by cluster","\n","\n",sep=""))
    dev.off()



    # Descriptive tables of the clusters
    clusterDesc=matrix(0,nrow=9,ncol=nbClust)
    for(i in 1:nbClust){
      clusterDesc[,i]=c(i,
                        length(which(cumsum(t(summaryClusters["Mean",,i])[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)])<50))+1,
                        length(which(cumsum(t(summaryClusters["Mean",,i])[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)])<90))+1,
                        length(which(t(summaryClusters["Median",,i])>50)),
                        length(which(resval[,i]>1.96)),
                        length(which(resval[,i]>3.29)),
                        length(which(apply(datSpecies,2,function (x) (clusters$clusinfo[i,1]-length(which(x[clusters$clustering==i]==0)))/clusters$clusinfo[i,1]*100)>50)),
                        length(which(apply(datSpecies,2,function (x) (clusters$clusinfo[i,1]-length(which(x[clusters$clustering==i]==0)))/clusters$clusinfo[i,1]*100)>90)),
                        clusters$clusinfo[i,1])
    }
    rownames(clusterDesc)=c("Number of species",
                            "to have 50% of catch", "to have 90% of catch",
                            "with a median higher than 50",
                            "with a test-value > 1.96", "with a test-value > 3.29",
                            "catch in 50% of the logevents", "catch in 90% of the logevents",
                            "Clusters size")
    colnames(clusterDesc)=1:nbClust
    clusterDesc2=as.data.frame(clusterDesc)


    # Summary tables of the clusters
    namesSpecies=matrix(NA,nrow=nbClust,ncol=10)
    namesCapt=matrix(NA,nrow=nbClust,ncol=5)
    nbSpeciesCatch = min(5,dim(t(summaryClusters["Mean",,]))[[2]])
    namesTarget=matrix(NA,nrow=nbClust,ncol=5)
    nbSpeciesVT = min(5,dim(target$tabnomespcib)[[2]])
    tabLibname=matrix(NA,nrow=nbClust,ncol=10)
    listLibname=list()

    for(i in 1:nbClust){
      namesCapt[i,]=colnames(t(summaryClusters["Mean",,i]))[order(t(summaryClusters["Mean",,i]),decreasing=TRUE)][1:nbSpeciesCatch]
      a=as.data.frame(t(summaryClusters["Mean",target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])],i]))
      colnames(a)= target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])]
      if(length(a)!=0){
        namesTarget[i,1:length(target$tabnomespcib[i,1:nbSpeciesVT][!is.na(target$tabnomespcib[i,1:nbSpeciesVT])])]=colnames(a[order(a,decreasing=TRUE)])
      }
      namesSpecies[i,1:length(union(namesCapt[i,],namesTarget[i,]))]=union(namesCapt[i,],namesTarget[i,])
    }
    
    for(i in 1:nbClust){
        listLibname[[i]]=lapply(as.list(namesSpecies[i,]), function(x) if(length(which(correspLevel7to5[,"X3A_CODE"]==x))==0) "NA" 
                                                    else correspLevel7to5[which(correspLevel7to5[,"X3A_CODE"]==x),"French_name"])
        tabLibname[i,]=unlist(lapply(listLibname[[i]], function(x) as.character(unlist(x))))
    }
         
    tabPropCatch=matrix(NA,nrow=nbClust,ncol=10)
    tabTestVal=matrix(NA,nrow=nbClust,ncol=10)
    tabPropLog=matrix(NA,nrow=nbClust,ncol=10)

    for(i in 1:nbClust){
      print("-----------------------------------------------------------------")
      print(paste("Cluster",i))
      propCatch=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)>=0.1)]
      tabPropCatch[i,1:length(propCatch)]=propCatch
      print(propCatch)
      testVal=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) resval[x,i]),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) t(summaryClusters["Mean",x,i])),digits=1)>=0.1)]
      tabTestVal[i,1:length(testVal)]=testVal
      print(testVal)
      propLog=round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) (clusters$clusinfo[i,1]-length(which(datSpecies[clusters$clustering==i,x]==0)))/clusters$clusinfo[i,1]*100),digits=1)[which(round(sapply(namesSpecies[i,][!is.na(namesSpecies[i,])],function(x) (clusters$clusinfo[i,1]-length(which(datSpecies[clusters$clustering==i,x]==0)))/clusters$clusinfo[i,1]*100),digits=1)>=0.1)]
      tabPropLog[i,1:length(propLog)]=propLog
      print(propLog)
    }

    tabClusters=array(0,dim=c(10,5,nbClust))
    dimnames(tabClusters)[[2]]=c("Libname","FAO","Test-value","% Catch","% Logevents")
    dimnames(tabClusters)[[3]]=paste("Cluster",1:nbClust)
    for(i in 1:nbClust){
      tabClusters[,,i]=cbind(tabLibname[i,],namesSpecies[i,],tabTestVal[i,],tabPropCatch[i,],tabPropLog[i,])
    }

    sizeTabClusters=numeric()
    for(i in 1:nbClust){
      sizeTabClusters[i]=min(length(namesSpecies[i,!is.na(namesSpecies[i,])]),length(tabPropCatch[i,!is.na(tabPropCatch[i,])]),length(tabTestVal[i,!is.na(tabTestVal[i,])]),length(tabPropLog[i,!is.na(tabPropLog[i,])]))
    }
    
    
    # Target Species
    # Intersection of species from tabClusters having : - % Cumulated Catch > thresholdCatch
    #                                                   - Test-value > thresholdTestValue
    #                                                   - % Logevents > thresholdLogevents 
    thresholdCatch=75
    thresholdTestValue=3
    thresholdLogevents=30
    
    sppCumCatch=list()
    sppTestValue=list()
    sppLogevents=list()
    targetSpeciesByCluster=list()
    
    for (i in 1:nbClust){
      percCatchCum=cumsum(as.numeric(tabClusters[,"% Catch",i]))
      nbSpSel=length(which(percCatchCum<thresholdCatch))+1
      sppCumCatch[[i]]=tabClusters[1:nbSpSel,"FAO",i]
      
      sppTestValue[[i]]=tabClusters[which(as.numeric(tabClusters[,"Test-value",i])>thresholdTestValue),"FAO",i]
      
      sppLogevents[[i]]=tabClusters[which(as.numeric(tabClusters[,"% Logevents",i])>thresholdLogevents),"FAO",i]
    
      targetSpeciesByCluster[[i]]=intersect(sppCumCatch[[i]],sppTestValue[[i]])
      targetSpeciesByCluster[[i]]=intersect(targetSpeciesByCluster[[i]],sppLogevents[[i]])
    }

    # List of metiers (level 7)
    listMetiersL7=list()
    for (i in 1:nbClust){
      metiersClusteri=targetSpeciesByCluster[[i]]
      metiersClusteri=as.character(unique(unlist(metiersClusteri)))
      metiersClusteri=paste(unlist(strsplit(metiersClusteri," ")),collapse=" ")
      listMetiersL7[[i]]=metiersClusteri
    }
  
    # Metier (level 7) of each logevent
    metierByLogeventL7=unlist(sapply(clusters$clustering,function(x) listMetiersL7[[x]]))



    # Create csv tables
    write.table(clusterDesc2,file="descClusters.csv",col.names=FALSE,sep=";")
    
    dfClust=data.frame()
    dfClust=paste("Clust ",1:nbClust,sep="")
    for(i in 1:nbClust){
      write.table(dfClust[i],file="tabClusters.csv",append=TRUE,row.names=FALSE,sep=";")
      tabClusti=as.data.frame(tabClusters[1:sizeTabClusters[i],,i])
      write.table(tabClusti,file="tabClusters.csv",append=TRUE,row.names=FALSE,sep=";")  
    }     


    LE_ID_clust=data.frame(LE_ID=LE_ID,clust=metierByLogeventL7)
    print(" --- end of step 3 ---")
    print(Sys.time()-t1)

    return(list(LE_ID_clust=LE_ID_clust, clusters=clusters,
    betweenVarClassifOnTot=betweenVarClassifOnTot, nbClust=nbClust,
    summaryClusters=summaryClusters, testValues=resval, 
    testValuesSpecies=target$tabnomespcib, percLogevents=percLogevents, 
    descClusters=clusterDesc2, tabClusters=tabClusters,
    targetSpecies=targetSpeciesByCluster))                                   

  }  else stop("methMetier must be hac, kmeans, pam or clara")
  # end of the methods


} # end of the function "getMetierClusters"
                          
