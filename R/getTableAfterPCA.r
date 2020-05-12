################################################################################
#  STEP 2 OF THE MULTIVARIATE CLASSIFICATION :                                 #
#         RUN A PCA ON THE DATASET FROM STEP 1                                 #
#         (2 CRITERIA : 70PERCENTS AND SCREETEST ARE AVAILABLE)                #
#         IT'S POSSIBLE TO KEEP THE DATASET FROM STEP 1 BY CHOOSING "NOPCA"    #
################################################################################


getTableAfterPCA = function(datSpecies,analysisName="",pcaYesNo="pca",criterion="70percents"){

  LE_ID <- rownames(datSpecies)
  NbSpecies <- dim(datSpecies)[2]
  datSpecies <- as.matrix(datSpecies,ncol=NbSpecies,nrow=length(LE_ID))

  print("######## STEP 2 PCA/NO PCA ON CATCH PROFILES ########")

  t1 <- Sys.time()
  print(paste("--- Selected method :",pcaYesNo, "---"))


  if(pcaYesNo=="pca"){
  
    print("Running PCA on all axes...")
    
    # PCA (Principal Component Analysis)
    log.pca = PCA(datSpecies, graph=FALSE, ncp=ncol(datSpecies))
    
    X11(5,5)
    plot.PCA(log.pca, choix = "var", axes = 1:2, new.plot=FALSE, title="", lim.cos2.var = 0.1)
    savePlot(filename=paste(analysisName,'species_projection_on_the_1_and_2_factorial_axis',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot.PCA(log.pca, choix = "var", axes = 1:2, new.plot=FALSE, title="", lim.cos2.var = 0.1)
    savePlot(filename=paste(analysisName,'species_projection_on_the_1_and_2_factorial_axis',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot.PCA(log.pca, choix = "var", axes = 2:3, new.plot=FALSE, title="", lim.cos2.var = 0.1)
    savePlot(filename=paste(analysisName,'species_projection_on_the_2_and_3_factorial_axis',sep="_"), type='png')
    dev.off()
    
    X11(5,5)
    plot.PCA(log.pca, choix = "ind", axes = 1:2, habillage = "ind", title="", new.plot=FALSE, cex=1.1)
    savePlot(filename=paste(analysisName,'projection_of_individuals_on_the_first_two_factorial_axis',sep="_"), type='png')
    dev.off()
 

    # Determine the number of axis to keep
    if(criterion=="70percents"){
      nbaxes=which(log.pca$eig[,3]>70)[1]   # we are taking the axis until having 70% of total inertia
      cat("--- number of axes:",nbaxes,"\n")
      cat("--- percentage inertia explained:",log.pca$eig[nbaxes,3],"\n")
    } else
    # OR
    if(criterion=="screetest"){
      nbaxes=which(scree(log.pca$eig[,1])[,"epsilon"]<0)[2]  # thanks to the scree-test
      cat("--- number of axes:",nbaxes,"\n")
      cat("--- percentage inertia explained:",log.pca$eig[nbaxes,3],"\n")
    } else stop("Criterion for PCA must be 70percents or screetest")


    # Eigenvalues and relative graphics
    log.pca$eig

    png(paste(analysisName,"Eigen values.png",sep="_"), width = 1200, height = 800)
    x=1:length(log.pca$eig[,1])
    barplot(log.pca$eig[,1],names.arg=x, main="Eigen values")
    dev.off()
   
    png(paste(analysisName,"Percentage of Inertia.png",sep="_"), width = 1200, height = 800)
    color=rep("grey",length(log.pca$eig[,1]))
    if(criterion=="screetest") color[1:nbaxes]="green"
    barplot(log.pca$eig[,2],names.arg=x, col=color, main="Percentage of Inertia of factorial axis", xlab="Axis", ylab="% of Inertia")
    dev.off()
   
    X11(5,5)
    op <- par(no.readonly = TRUE) 
    par(mar=c(4,4,1,1))
    color=rep("grey",length(log.pca$eig[,1]))
    if(criterion=="70percents") color[1:nbaxes]="green"
    barplot(log.pca$eig[,3], col=color, ylab="", xlab="", axes=FALSE, cex.names=2)
    axis(2,  las=2)
    abline(h=70, col="red")
    text(1,72, "70% of Inertia", col = "red", adj = c(0, -.1))
    mtext("Axes", side=1, adj=0.5, outer=FALSE, line=+1, font=1, cex=1.5)
    mtext("% of Inertia", side=2, adj=0.5, outer=FALSE, line=+2.5, font=1, cex=1.5)
    savePlot(filename = paste(analysisName,"Cumulative Percentage of Inertia.png",sep="_"),type ="png")
    par(op)
    dev.off()
   

    # Projection of variables "species" on the first factorial axis
    png(paste(analysisName,"Projection of Species on first factorial axis.png",sep="_"), width = 1200, height = 800)
    op <- par(mfrow=c(2,3))
    plot(log.pca,choix="var",axes = c(1, 2),new.plot=FALSE,lim.cos2.var = 0.3)
    plot(log.pca,choix="var",axes = c(2, 3),new.plot=FALSE,lim.cos2.var = 0.3)
    plot(log.pca,choix="var",axes = c(1, 3),new.plot=FALSE,lim.cos2.var = 0.3)
    plot(log.pca,choix="var",axes = c(1, 4),new.plot=FALSE,lim.cos2.var = 0.3)
    plot(log.pca,choix="var",axes = c(2, 4),new.plot=FALSE,lim.cos2.var = 0.3)
    #plot(log.pca,choix="var",axes = c(3, 4),new.plot=FALSE,lim.cos2.var = 0.3)
    par(op)
    title(main=paste("Projection of Species on first factorial axis","\n","\n",sep=""))
    dev.off()


    # PCA with the good number of axis
    print("Retaining Principal Components of selected axes...")
    log.pca=log.pca$ind$coord[,1:nbaxes]
    datLog=round(log.pca,4)


  } else


  if(pcaYesNo=="nopca"){
    datLog=datSpecies
  }  else stop("pcaYesNo must be pca or nopca")

  Store(objects())
  gc(reset=TRUE)

  print(" --- end of step 2 ---")
  print(Sys.time()-t1)

  return(datLog)

}

