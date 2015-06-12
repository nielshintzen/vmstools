# Assign one métier to each logbook event.

# Authors #

Nicolas Deporte, Sébastien Demanèche, Stéphanie Mahévas (IFREMER, France)

Clara Ulrich, François Bastardie (DTU Aqua, Denmark)



# Introduction #

The goal of this part was to investigate a range of multivariate statistical methods and to develop a tool selecting the most appropriate one for assigning one métier to each logbook event. This tool allows defining métiers at DCF Level7 and aggregates them into the DCF Level5.

Logbooks data are the main source of data for getting information on the fishing activities. From logbooks data, information of EU level gear type and selectivity feature of the gear are directly available. However, fishers do not declare which species they are actually targeting when fishing, and therefore the DCF Level5 and the more precise DCF Level7, expected to reflect the fishing intention, has to be inferred from the available data, and in particular the catch data. Following, the ICES (2003) recommendation, analyses have been performed on landings composition expressed in cash value as this may reflect more accurately the actual targeting choices of the fishermen.

The work consisted in designing several sequential steps to have a complete tool firstly allowing identifying métiers from logbooks raw data and then assigning a métier to a new dataset.
These sequential steps aim at:

i) Identifying the main species out of all species recorded and reducing the dataset to these key species only (Step 1). These key species are selected from the application of three different filtering methods.

ii) Running or not a PCA to increase the quality of the final classification (PCA is likely to discard potential irrelevant information) (Step 2).

iii) Running a selection of clustering methods (HAC or K-Means or CLARA algorithm) to group similar logbook events and defining the characteristic species-based assemblage of each cluster (DCF Level7) (Step 3).

iv) Converting them to a corresponding DCF Level5 category and comparing with simpler ordination methods (Step 4).

v) Returning a complete Eflalo dataset including the métier linked to each logbook event (Step 5).

vi) Predicting the classification of any new logbook event into the categories of defined métiers (Step 6).



# Functions #

This analysis needs the following functions to run:

  * multivariateClassificationUtils.r
  * selectMainSpecies.r
  * extractTableMainSpecies.r
  * getTableAfterPCA.r
  * getMetierClusters.r
  * compareToOrdination.r
  * predictMetier.r
  * getEflaloMetierLevel7.r




# Details #

```

  data(eflalo)
  data(correspLevel7to5)
  data(correspMixedMetier)

  # Initialization
  # Choose the working directory
  path="C:/metier_analysis"
  analysisName="metier_analysis_OTB"  

  # Settings for the step 1
  runHACinSpeciesSelection=TRUE
  paramTotal=95
  paramLogevent=100
  # Setting for the step 2
  critPca="PCA_70"
  # Setting for the step 3
  algoClust="CLARA"

  # Format
  eflalo=formatEflalo(eflalo)
  eflalo=eflalo[eflalo$LE_GEAR=="OTB",]

  #-------------------------------------------------------------------------
  # I. GETTING THE DATA IN AND CLEANING FOR MISSING AND NEGATIVE VALUES ETC
  #-------------------------------------------------------------------------
  
  print("--- CREATING DIRECTORIES AND REDUCING THE EFLALO DATASET TO THE ONLY DATA NECESSARY FOR THE ANALYSIS ---")
  cat("\n") 
    
  # Creating the directory of the analysis
  if (!file.exists(path)) dir.create(path)
  setwd(path)
  if (!file.exists(analysisName)) dir.create(analysisName)
  setwd(file.path(path,analysisName))
  # Delete old R cache
  if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)                                            
  eflalo_ori=eflalo # Keeping this in cached memory for making the final merging at the end
  Store(eflalo_ori)

  # ! KEEPING ONLY LE_ID AND THE OUTPUT YOU WANT TO GET  (KG/EURO)
  eflalo=eflalo[,c("LE_ID",grep("EURO",names(eflalo),value=T))]
  eflalo[is.na(eflalo)]=0

  # Removing negative and null values
  null.value=vector()
  for (i in grep("EURO",names(eflalo))) null.value=c(null.value,which(eflalo[,i]<0))
  null.value=c(null.value,which(apply(eflalo[,2:ncol(eflalo)],1,sum,na.rm=T)==0))

  if(length(null.value)!=0) {LogEvent.removed=eflalo[sort(unique(null.value)),] ; eflalo=eflalo[-sort(unique(null.value)),]}

  # Rename species names
  names(eflalo)[-1]=unlist(lapply(strsplit(names(eflalo[,-1]),"_"),function(x) x[[3]]))

  # Removing miscellaneous species
  eflalo=eflalo[,!names(eflalo)=="MZZ"]
  save(eflalo, file="eflalo_cleaned.Rdata")


  #----------------------------------------------------------------------------------------------------------
  # II. STEP 1 - EXPLORING THE VARIOUS METHODS FOR IDENTIFYING MAIN SPECIES AND KEEPING THEM IN THE DATA SET
  #----------------------------------------------------------------------------------------------------------
  
  print("--- EXPLORING THE DATA FOR SELECTION OF MAIN SPECIES ---")
  cat("\n") 

  # Exploration of main species
  explo=selectMainSpecies(eflalo,analysisName,RunHAC=runHACinSpeciesSelection,DiagFlag=FALSE)

  # Step 1 : selection of main species
  Step1=extractTableMainSpecies(eflalo,explo$namesMainSpeciesHAC,paramTotal=paramTotal,paramLogevent=paramLogevent)

  save(explo,Step1,file="Explo_Step1.Rdata")
  
  
  #-----------------------------
  # III. STEP 2 - PCA - NO-PCA
  #-----------------------------

  if (!file.exists(critPca)) dir.create(critPca)
  setwd(file.path(path,analysisName,critPca))
  if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)

  rowNamesSave=row.names(Step1)
  row.names(Step1)=1:nrow(Step1)

  if (critPca=="PCA_70") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="70percents") else    # criterion="70percents"
  if (critPca=="PCA_SC") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="screetest") else    # criterion="screetest"
  if (critPca=="NO_PCA") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="nopca",criterion=NULL)

  row.names(Step1)=rowNamesSave
  row.names(Step2)=rowNamesSave

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
  
  compMetiers=compareToOrdination(dat=eflalo,Step2=Step2,clusters=clusters,targetSpecies=Step3$targetSpecies)
  
  save(compMetiers,file="compMetiers.Rdata")
  
  
  #-------------------------------------
  # VI. STEP 5 - MERGING BACK TO EFLALO
  #-------------------------------------

  # Choosing the final option
  setwd(file.path(path,analysisName))

  if(!nrow(eflalo)==nrow(Step3$LE_ID_clust)) print("--error : number of lines in step 3 not equal to input eflalo, please check!!--")

  eflalo=cbind(eflalo,CLUSTER=Step3$LE_ID_clust[,"clust"])

  # Now reload the full data set
  eflalo_ori[-sort(unique(null.value)),"CLUSTER"]=Step3$LE_ID_clust[,"clust"]

  print("Congratulation !! You have now a fully working eflalo dataset with a metier Level 7 !")
  

  #-------------------------------------------------------------------------------------------------------------
  # VII. STEP 6 - PREDICTING THE CLASSIFICATION OF ANY NEW LOGBOOK EVENT INTO THE CATEGORIES OF DEFINED MÉTIERS
  #-------------------------------------------------------------------------------------------------------------

  data(eflalo)

  # Format
  newEflalo=formatEflalo(eflalo)
  newEflalo=newEflalo[newEflalo$LE_GEAR=="OTB",]

  # GETTING THE NEW DATA IN AND CLEANING FOR MISSING AND NEGATIVE VALUES ETC 
  datPred=newEflalo[,c("LE_ID",grep("EURO",names(newEflalo),value=T))]
  rm(newEflalo)
  gc()
  datPred[is.na(datPred)]=0
    
  # Removing negative and null values
  null.value=vector()
  for (i in grep("EURO",names(datPred))) null.value=c(null.value,which(datPred[,i]<0))
  null.value=c(null.value,which(apply(datPred[,2:ncol(datPred)],1,sum,na.rm=T)==0))
    
  if(length(null.value)!=0) {LogEvent.removed=datPred[sort(unique(null.value)),] ; datPred=datPred[-sort(unique(null.value)),]} 
  
  # Rename species names
  names(datPred)[-1]=unlist(lapply(strsplit(names(datPred[,-1]),"_"),function(x) x[[3]]))
    
  # Removing miscellaneous species
  datPred=datPred[,!names(datPred)=="MZZ"]

  # Create a new folder for Step 6
  option_step6="predict_Metiers"
  setwd(file.path(path,analysisName,critPca,algoClust))
  if (!file.exists(option_step6)) dir.create(option_step6)
  setwd(file.path(path,analysisName,critPca,algoClust,option_step6))
  if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)

  # Prepare datPred
  le_id_datPred=datPred[,"LE_ID"]
  datPred=datPred[,-1]
  datPred=as.matrix(datPred)
  rownames(datPred)=le_id_datPred

  if (algoClust=="HAC")     clustersAffectation=Step3$clusters
  if (algoClust=="CLARA")   clustersAffectation=Step3$clusters$clustering
  if (algoClust=="KMEANS")  clustersAffectation=Step3$clusters$cluster  

  # Predict the metiers of the new logbook events of newData from the metiers found in Step 3
  metierPred=predictMetier(learningData=Step1,clustersAffectation=clustersAffectation,newData=datPred)
  
  save(metierPred,file="metierPred.Rdata")

```



**Note** : The steps 1 to 5 (included) can be summarized by the function getEflaloMetierLevel7. This latter returns a fully working eflalo dataset with a metier at DCF Level7 for each logbook event, accompanied by several graphics and results about the classification.


```

  data(eflalo)
 
  # Choose the working directory
  path="C:/metier_analysis2"
  analysisName="metier_analysis_OTB"  
  if (!file.exists(path)) dir.create(path)
  setwd(path)
  if (!file.exists(analysisName)) dir.create(analysisName)
  setwd(file.path(path,analysisName))

  # Format
  eflalo=formatEflalo(eflalo)
  eflalo=eflalo[eflalo$LE_GEAR=="OTB",]

  # Return a fully working eflalo dataset with a metier at DCF Level7 for each logbook event
  eflalo_metiers_level7=getEflaloMetierLevel7(eflalo, analysisName, path, critData="EURO", runHACinSpeciesSelection=TRUE, paramTotal=95, paramLogevent=100, critPca="PCA_70", algoClust="CLARA")

  save(eflalo_metiers_level7,file="eflalo_metiers_level7.Rdata")

```