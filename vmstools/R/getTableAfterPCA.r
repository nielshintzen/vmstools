################################################################################
#  STEP 2 OF THE MULTIVARIATE CLASSIFICATION :                                 #
#         RUN A PCA ON THE DATASET FROM STEP 1                                 #
#         (2 CRITERIA : 70PERCENTS AND SCREETEST ARE AVAILABLE)                #
#         IT'S POSSIBLE TO KEEP THE DATASET FROM STEP 1 BY CHOOSING "NOPCA"    #
################################################################################




#' Finding Metiers from a reduced EFLALO dataset, step 2: Options for running a
#' PCA on the selected species.
#' 
#' This function represents the second step in the multivariate analysis of
#' logbooks data for identifying metiers.
#' 
#' This step of the classification is a user-driven choice of running or not a
#' Principal Component Analysis (PCA) on the catch per logevents of the species
#' retained, expressed in percentage. A PCA will help reducing the
#' multi-dimensional catch matrix to a smaller number of informative components
#' represented by the first n- axes of the PCA transformation. The number of
#' axes can be chosen either by using a scree test (criterion 'screetest')
#' looking for the significant marginal increases of explained inertia, or by
#' selecting all axes cumulating 70 percent of explained inertia (criterion
#' '70percents').
#' 
#' 
#' @param datSpecies numerical matrix with Logevents as lines and species as
#' columns, with percentage values (between 0 and 100) of each species in the
#' logevent catches.  Logevent ID (LE_ID) should be as row names. Typically,
#' this table will be produced from a eflalo dataset using the function
#' extractTableMainSpecies()
#' @param analysisName character, the name of the run. Used for the file name
#' of the plots.
#' @param pcaYesNo character. An indication of whether a PCA should be
#' performed. Use "pca" if a PCA should be run, and "nopca" if it shouldn't.
#' @param criterion character. Criterion used for selected the number of
#' principal components (axes) retained.  The number of axes can be chosen
#' either by using a scree test (criterion 'screetest'), looking for the
#' second-order unsignificant marginal increases of explained inertia, or by
#' selecting all axes cumulating up to 70 percent of explained inertia
#' (criterion '70percents'). If pcaYesNo="nopca", then this criterion should be
#' set as NULL. Experience has shown that the '70percents' criterion may often
#' select a larger number axes than the scree test criterion, and may therefore
#' be more appropriate for analysing large and heterogeneous datasets, by
#' retaining more information for the subsequent clustering.
#' @return If a PCA is run, a number of graphs are produced and directly saved
#' in the working directory. These graphs describe 1) species projection on the
#' first factorial axes, 2) individuals (logevents) projections on the first
#' factorial axes, 3) actual and cumulative percentage of inertia of factorial
#' axes, and 4) eigenvalues.
#' 
#' If a PCA is run, the function returns a numerical matrix with Logevents as
#' lines and selected Principal Components as columns.  If no PCA is run, the
#' function returns the same matrix as the input, with percentage values by
#' species.
#' @note A number of libraries are initially called for the whole metier
#' analyses and must be installed :
#' (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)
#' @author Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER,
#' France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)
#' @seealso selectMainSpecies(), extractTableMainSpecies()
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
#'                             criterion="70percents")
#' 
#'   row.names(Step1) <- rowNamesSave
#'   row.names(Step2) <- rowNamesSave
#'   
#'   }
#' 
#' 
#' @export getTableAfterPCA
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

