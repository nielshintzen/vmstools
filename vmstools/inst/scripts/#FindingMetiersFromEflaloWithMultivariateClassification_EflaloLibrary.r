###################################################################################
##     FINDING METIERS FROM EFLALO USING MULTIVARIATE CLUSTERING METHODS         ##
##                                                                               ##
##              LOT 2 - VMS LOGBOOKS (WP2)                                       ##
##                                                                               ##
##                                                                               ##
## Authors : Nicolas Deporte (IFREMER / OCEANIC DEVELOPMENT)                     ##
##           Stephanie Mahévas, Sebastien Demaneche (IFREMER)                    ##
##           Clara Ulrich, Francois Bastardie (DTU Aqua)                         ##
##                                                                               ##
## Last update : October 2011,                                                   ##
##                                                                               ##
## Runs with R 2.13.2                                                            ##
##                                                                               ##
###################################################################################



rm(list=ls(all=TRUE))
gc(reset=TRUE)

path = "Q:\\Nicolas\\Scripts\\R\\Analyses\\Test"
setwd(path) # you must choose the path of your working directory

memory.limit(size=4000)

#-----------------------------
# I. GETTING THE DATA IN AND CLEANING FOR MISSING AND NEGATIVE VALUES ETC
#-----------------------------
 
country = "All"
year = 2007
AreaCodename = "3a4"
Gear = c("OTB")

analysisName=paste(country,"_",Gear,year,"_",AreaCodename,sep="")

data(eflalo) # load a data eflalo dataset e.g. the one given as illustration in vmstools
#load("eflalo.rda")
dat = eflalo

# creating the directory of the analysis
if (!file.exists(analysisName)) dir.create(analysisName)
setwd(paste(path,analysisName,sep="/"))
#delete old R cache
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)      


eflalo_ori = dat # keeping this in cached memory for making the final merging at the end
Store(eflalo_ori)

# ! KEEPING ONLY LE_ID AND THE OUTPUT YOU WANT TO GET  (KG/EURO)
dat = dat[,c("LE_ID",grep("EURO",names(dat),value=T))]
dat[is.na(dat)]=0

#removing negative and null values
null.value = vector()
for (i in grep("EURO",names(dat))) null.value = c(null.value,which(dat[,i]<0))
null.value = c(null.value,which(apply(dat[,2:ncol(dat)],1,sum,na.rm=T)==0))

if(length(null.value)!=0) {
  LogEvent.removed = dat[sort(unique(null.value)),]
  dat = dat[-sort(unique(null.value)),]
}

names(dat)[-1]=unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))

#removing miscellaneous species
dat = dat[,!names(dat)=="MZZ"]

save(dat, file="dat_2007.Rdata")
#load("dat_2007.Rdata")

#-----------------------------
# II. EXPLORING THE VARIOUS METHODS FOR IDENTIFYING MAIN SPECIES AND KEEPING THEM IN THE DATA SET (STEP 1)
#-----------------------------
#EXPLORATION
explo=selectMainSpecies(dat,analysisName,RunHAC=TRUE,DiagFlag=FALSE)

# Step 1 : selection of main species
Step1=extractTableMainSpecies(dat,explo$namesMainSpeciesHAC,paramTotal=95,paramLogevent=100)

save(explo,Step1,file="Explo_Step1.Rdata")

#load("Explo_Step1.Rdata")

#-----------------------------
# III. STEP 2 - PCA - NO-PCA
#-----------------------------

# Step 2 : PCA

rownames(Step1)=paste("eflalo",1:nrow(Step1),sep="")     # to fix the problem : several lines have the same name

option_step2="PCA_70"

setwd(paste(path,analysisName,sep="/"))
if (!file.exists(option_step2)) dir.create(option_step2)
setwd(paste(path,analysisName,option_step2,sep="/"))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)      


if (option_step2=="PCA_70") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="70percents") else    # criterion="70percents"
if (option_step2=="PCA_SC") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="screetest") else    # criterion="screetest"
if (option_step2=="NO_PCA") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="nopca",criterion=NULL)    

save(Step2,file="Step2.Rdata")
                                      
#load("Step2.Rdata")

#-----------------------------
# IV. STEP 3 - CLUSTERING METHOD : HAC, CLARA OR KMEANS
#-----------------------------
data(correspLevel7to5)
data(correspMixedMetier)

option_step3="CLARA"

setwd(paste(path,analysisName,option_step2,sep="/"))
if (!file.exists(option_step3)) dir.create(option_step3)
setwd(paste(path,analysisName,option_step2,option_step3,sep="/"))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)      


if (option_step3=="HAC")    Step3=getMetierClusters(Step1,Step2,analysisName,methMetier="hac",param1="euclidean",param2="ward") else    
if (option_step3=="CLARA")  Step3=getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="clara",param1="euclidean",param2=NULL) else    
if (option_step3=="KMEANS") Step3=getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="kmeans",param1=NULL,param2=NULL)    

save(Step3,file="Step3.Rdata")
#load("Step3.Rdata")


#-----------------------------
# V. STEP 4 - MERGING BACK TO EFLALO
#-----------------------------

##choosing the final option
#setwd(paste(path,analysisName,sep="/"))
#
#load("PCA_70/CLARA/Step3.Rdata")

if(!nrow(dat)==nrow(Step3$LE_ID_clust)) print("--error : number of lines in step 3 not equal to input eflalo, please check!!--")

dat_WithClusters = cbind(dat,CLUSTER=Step3$LE_ID_clust[,"clust"])

#now reload the full data set
if(length(null.value)!=0){ eflalo_ori[-sort(unique(null.value)),"CLUSTER"] = Step3$LE_ID_clust[,"clust"]
} else { eflalo_ori[,"CLUSTER"] = Step3$LE_ID_clust[,"clust"] }
eflalo_ori_WithClusters=eflalo_ori

save(dat_WithClusters,file="dat_WithClusters.Rdata")
save(eflalo_ori_WithClusters,file="eflalo_ori_WithClusters.Rdata")


#-----------------------------------------------------
# VI. STEP 5 - COMPARISON WITH ORDINATION METHODS
#-----------------------------------------------------

# load previous R objects (Step1,Step2,Step3)
#setwd(paste(path,analysisName,sep="/"))
#load("dat_2007.Rdata")
#load("Explo_Step1.Rdata")
#option_step2="PCA_70"
#setwd(paste(path,analysisName,option_step2,sep="/"))
#load("Step2.Rdata")
#option_step3="CLARA"
#setwd(paste(path,analysisName,option_step2,option_step3,sep="/"))
#load("Step3.Rdata")
compOrdin="CompOrdin"
if (!file.exists(compOrdin)) dir.create(compOrdin)
setwd(paste(path,analysisName,option_step2,option_step3,compOrdin,sep="/"))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)
compMetiers=compareToOrdination(dat=dat,Step2=Step2,clusters=Step3$clusters$clustering,targetSpecies=Step3$targetSpecies)
save(compMetiers,file="compMetiers.Rdata")


#-------------------------------------------------------------------------------------------
# VII. STEP 6 - Predicting Metier for current year using clustering performed previous year
#-------------------------------------------------------------------------------------------

# load previous R objects (Step1, Step3)
#setwd(path)
#getwd()
#option_step2="PCA_70"
#option_step3="CLARA"
#load(paste(analysisName,"/",option_step2,"/",option_step3,"/Step3.Rdata",sep=""))
clust2007=Step3$clusters$clustering
#rm(Step3)
#gc()
#load(paste(analysisName,"/Explo_Step1.Rdata",sep=""))
#rm(explo)
#gc()

# load your new dataset (called dat here)
#! KEEPING ONLY LE_ID AND THE OUTPUT YOU WANT TO GET  (KG/EURO)
new_country = "All"
new_year = 2008
new_AreaCodename = "3a4"
new_Gear = c("OTB")
newAnalysisName=paste(new_country,"_",new_Gear,new_year,"_",new_AreaCodename,sep="")

#load("dat_2008.Rdata")
dat = eflalo

## if dat is ready, use :
#    datPred=dat
#    rm(dat)
#    gc()

# else, if dat is a EFLALO format dataset, use :
    datPred = dat[,c("LE_ID",grep("EURO",names(dat),value=T))]
    rm(dat)
    gc()
    datPred[is.na(datPred)]=0
    
    #removing negative and null values
    null.value = vector()
    for (i in grep("EURO",names(datPred))) null.value = c(null.value,which(datPred[,i]<0))
    null.value = c(null.value,which(apply(datPred[,2:ncol(datPred)],1,sum,na.rm=T)==0))
    
    if(length(null.value)!=0) {LogEvent.removed = datPred[sort(unique(null.value)),] ; datPred = datPred[-sort(unique(null.value)),]}
    #Store(LogEvent.removed)
    
    names(datPred)[-1]=unlist(lapply(strsplit(names(datPred[,-1]),"_"),function(x) x[[3]]))
    
    #removing miscellaneous species
    datPred = datPred[,!names(datPred)=="MZZ"]


# and continue here :
#create a new folder for step6
option_step5 = "Predict_2008"
setwd(paste(path,analysisName,option_step2,option_step3,sep="/"))
if (!file.exists(option_step5)) dir.create(option_step5)
setwd(paste(path,analysisName,option_step2,option_step3,option_step5,sep="/"))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)

# prepare datPred
le_id_datPred=datPred[,"LE_ID"]
datPred=datPred[,-1]
datPred=as.matrix(datPred)
rownames(datPred)=le_id_datPred

metierPred2008 = predictMetier(learningData=Step1,clustersAffectation=clust2007,newData=datPred)
save(metierPred2008,file="metierPred2008.Rdata")
