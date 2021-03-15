#' Finding Metiers from a full EFLALO dataset, general function calling all the
#' various steps of the analyses at once
#' 
#' This function will run the various steps of the multivariate analyses in the
#' right order, allowing the full anlysis to be run at once.
#' 
#' 
#' @param dat a data.frame in eflalo format. It should be previously subsetted
#' for the desired region and gear of interest.
#' @param analysisName character, the name of the run. Used for the file name
#' of the plots.
#' @param path character, the path of the main directory where all data and
#' plots will be stored. A number of sub-directories will be created by the
#' function.
#' @param critData character. Choice if the analysis must be performed on
#' eflalo data in weight ("KG") or in value ("EURO")
#' @param runHACinSpeciesSelection Boolean. the HAC can be time consuming for
#' large datasets, and can therefore be switched off for quicker runs.
#' @param paramTotal numeric. The percentage threshold for the Total method.
#' All species ranked in terms of percentage in total catches and cumulating up
#' to this percentage value will be retained. Default value is 95\%.
#' @param paramLogevent numeric. The percentage threshold for the Logevent
#' method. All species representing at least this threshold of at least one
#' logevent (one line) will be retained. Default value is 100\% (logevents
#' catching one and only one species).
#' @param critPca character. Choice if a PCA must be run or not (="NO_PCA"),
#' and if yes, which criterion is used for selecting the number of axes, to be
#' chosen between scree test (="PCA_SC") or 70\% of cumulated inertia
#' (="PCA_70")
#' @param algoClust character. Choice of the clustering algorithm, either
#' "CLARA", "HAC" or "KMEANS".
#' @return
#' 
#' The function returns the same eflalo data set as was input, with an
#' additional column giving the cluster name where each logevent has been
#' allocated.
#' 
#' Beside, a number of sub-directories are created, with intermediate results
#' and plots.
#' @note A number of libraries are initially called for the whole metier
#' analyses and must be installed :
#' (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)
#' @author Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER,
#' France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)
#' @seealso \code{\link{selectMainSpecies}},
#' \code{\link{extractTableMainSpecies}}, \code{\link{getTableAfterPCA}},
#' \code{\link{getMetierClusters}},
#' @references Development of tools for logbook and VMS data analysis. Studies
#' for carrying out the common fisheries policy No MARE/2008/10 Lot 2
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' data(eflalo)
#' 
#' eflalo <- formatEflalo(eflalo)
#' 
#' eflalo <- eflalo[eflalo$LE_GEAR=="OTB",]
#' 
#' # Note that output plots will be sent to getwd()
#' path <- getwd()
#' analysisName <- "metier_analysis_OTB"
#' 
#' # Return a fully working eflalo dataset with a metier at DCF Level7 for each
#' # logbook event
#' eflalo_metiers_level7 <- getEflaloMetierLevel7(eflalo, analysisName, path,
#'     critData="EURO",runHACinSpeciesSelection=TRUE, paramTotal=95,
#'     paramLogevent=100,critPca="PCA_70", algoClust="CLARA")
#'         
#' }
#' 
#' 
#' @export getEflaloMetierLevel7
getEflaloMetierLevel7=function(dat, analysisName, path, critData="EURO", runHACinSpeciesSelection=TRUE, paramTotal=95, paramLogevent=100, critPca="PCA_70", algoClust="CLARA"){

  #-------------------------------------------------------------------------
  # I. GETTING THE DATA IN AND CLEANING FOR MISSING AND NEGATIVE VALUES ETC
  #-------------------------------------------------------------------------

  # Load the table linking 3A-CODE (FAO CODE of species) to the species assemblage (level 5).
  data(correspLevel7to5)
  
  # Load the table linking mixed metiers (composed by 2 simple metiers) to their official code of mixed metiers level 5 (FAO - 3 characters).
  data(correspMixedMetier)
  
  timeStart=Sys.time()
  
  print("--- CREATING DIRECTORIES AND REDUCING THE EFLALO DATASET TO THE ONLY DATA NECESSARY FOR THE ANALYSIS ---")
  cat("\n") 
    
  # Creating the directory of the analysis
  if (!file.exists(analysisName)) dir.create(analysisName)
  setwd(file.path(path,analysisName))
  # Delete old R cache
  if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)                                            

  eflalo_ori = dat # Keeping this in cached memory for making the final merging at the end
  Store(eflalo_ori)

  # ! KEEPING ONLY LE_ID AND THE OUTPUT YOU WANT TO GET  (KG/EURO)
  dat=dat[,c("LE_ID",grep(critData,names(dat),value=TRUE))]
  dat[is.na(dat)]=0

  # Removing negative and null values
  null.value <- vector()
  for (i in grep(critData,names(dat))) null.value <- c(null.value,which(dat[,i]<0))
  null.value <- c(null.value,which(apply(dat[,2:ncol(dat)],1,sum,na.rm=TRUE)==0))

  if(length(null.value)!=0) {LogEvent.removed <- dat[sort(unique(null.value)),] ; dat <- dat[-sort(unique(null.value)),]}

  # Rename species names
  names(dat)[-1]=unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))

  # Removing miscellaneous species
  dat <- dat[,!names(dat)=="MZZ"]
  save(dat, file="dat_cleaned.Rdata")


  #----------------------------------------------------------------------------------------------------------
  # II. EXPLORING THE VARIOUS METHODS FOR IDENTIFYING MAIN SPECIES AND KEEPING THEM IN THE DATA SET (STEP 1)
  #----------------------------------------------------------------------------------------------------------
  
  print("--- EXPLORING THE DATA FOR SELECTION OF MAIN SPECIES ---")
  cat("\n") 

  # Exploration of main species
  explo=selectMainSpecies(dat,analysisName,RunHAC=runHACinSpeciesSelection,DiagFlag=FALSE)

  # Step 1 : selection of main species
  Step1=extractTableMainSpecies(dat,explo$namesMainSpeciesHAC,paramTotal=paramTotal,paramLogevent=paramLogevent)

  save(explo,Step1,file="Explo_Step1.Rdata")
  
  rowNamesSave <- row.names(Step1)
  row.names(Step1) <- 1:nrow(Step1)
  
  #-----------------------------
  # III. STEP 2 - PCA - NO-PCA
  #-----------------------------

  if (!file.exists(critPca)) dir.create(critPca)
  setwd(file.path(path,analysisName,critPca))
  if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)

  if (critPca=="PCA_70") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="70percents") else    # criterion="70percents"
  if (critPca=="PCA_SC") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="screetest") else    # criterion="screetest"
  if (critPca=="NO_PCA") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="nopca",criterion=NULL)

  row.names(Step1) <- rowNamesSave
  row.names(Step2) <- rowNamesSave
  
  save(Step1,file="Step1.Rdata")
  save(Step2,file="Step2.Rdata")
  
  #-------------------------------------------------------
  # IV. STEP 3 - CLUSTERING METHOD : HAC, CLARA OR KMEANS
  #-------------------------------------------------------

  if (!file.exists(algoClust)) dir.create(algoClust)
  setwd(file.path(path,analysisName,critPca,algoClust))
  if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)

  if (algoClust=="HAC")    Step3=getMetierClusters(Step1,Step2,analysisName,methMetier="hac",param1="euclidean",param2="ward") else
  if (algoClust=="CLARA")  Step3=getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="clara",param1="euclidean",param2=NULL) else
  if (algoClust=="KMEANS") Step3=getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="kmeans",param1=NULL,param2=NULL)

  save(Step3,file="Step3.Rdata")


  #------------------------------------------------
  # V. STEP 4 - COMPARISON WITH ORDINATION METHODS
  #------------------------------------------------

  compOrdin="CompOrdin"
  if (!file.exists(compOrdin)) dir.create(compOrdin)
  setwd(file.path(path,analysisName,critPca,algoClust,compOrdin))
  if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)
  
  if (algoClust=="HAC")     clusters=Step3$clusters
  if (algoClust=="CLARA")   clusters=Step3$clusters$clustering
  if (algoClust=="KMEANS")  clusters=Step3$clusters$cluster  
  
  compMetiers=compareToOrdination(dat=dat,Step2=Step2,clusters=clusters,targetSpecies=Step3$targetSpecies)
  save(compMetiers,file="compMetiers.Rdata")
  
  
  #-------------------------------------
  # VI. STEP 5 - MERGING BACK TO EFLALO
  #-------------------------------------

  # Choosing the final option
  setwd(file.path(path,analysisName))

  #load(paste(path,analysisName,critPca,algoClust,"Step3.Rdata",sep="/"))

  if(!nrow(dat)==nrow(Step3$LE_ID_clust)) print("--error : number of lines in step 3 not equal to input eflalo, please check!!--")

  dat <- cbind(dat,CLUSTER=Step3$LE_ID_clust[,"clust"])

  # Now reload the full data set
  if(length(null.value)==0){
    eflalo_ori[,"CLUSTER"]=Step3$LE_ID_clust[,"clust"]
  }else{
    eflalo_ori[-sort(unique(null.value)),"CLUSTER"]=Step3$LE_ID_clust[,"clust"]
  }
  
  cat("\n")
  print("Congratulation !! You have now a fully working eflalo dataset with a metier Level 7 !")
  cat("\n")
  print(Sys.time()-timeStart)
  
  return(eflalo_ori)
}
