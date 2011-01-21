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
## Last update : December 2010   ,                                               ##
##                                                                               ##
## Runs with R 2.11.1                                                            ##
##                                                                               ##
###################################################################################



rm(list=ls(all=TRUE))
gc(reset=TRUE)

#path <-"G:/DossiersASynchroniser/Stat-Methode-Mequapro/EtudePCP/R/RDec2010/"
path <- "C:/Nicolas/Scripts/R/Analyses"
setwd(path) # you must choose the path of your working directory

#source("programs/multivariateClassificationUtils.r")
#source("extractTableMainSpecies.r")
#source("programs/getTableAfterPCA.r")
#source("programs/getMetierClusters.r")
#source("programs/selectMainSpecies.r")
#source("programs/predictMetier.r")
#memory.limit(size=4000)

source("multivariateClassificationUtils.r")
source("extractTableMainSpecies.r")
source("getTableAfterPCA.r")
source("getMetierClusters.r")
source("selectMainSpecies.r")
source("compareMultivariateClassificationToOrdination.r")
source("level7to5.r")
source("level5.r")
source("predictMetier.r")
memory.limit(size=4000)



#-----------------------------
# I. GETTING THE DATA IN AND CLEANING FOR MISSING AND NEGATIVE VALUES ETC
#-----------------------------
# 
country <- "All"
year <- 2007
AreaCodename <- "3a4"
Gear <- c("OTB")

analysisName=paste(country,"_",Gear,year,"_",AreaCodename,sep="")

# load your own dataset (called dat1 here)
#load(paste("data/All_eflalo_2007OTB3a4.Rdata",sep=""))
load("All_eflalo_2007OTB3a4.Rdata")

# creating the directory of the analysis
if (!file.exists(analysisName)) dir.create(analysisName)
setwd(paste(path,analysisName,sep="/"))
#delete old R cache
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)      


eflalo_ori <- dat1 # keeping this in cached memory for making the final merging at the end
Store(eflalo_ori)

# ! KEEPING ONLY LE_ID AND THE OUTPUT YOU WANT TO GET  (KG/EURO)
dat1 <- dat1[,c("LE_ID",grep("EURO",names(dat1),value=T))]
dat1[is.na(dat1)]=0

#removing negative and null values
null.value <- vector()
for (i in grep("EURO",names(dat1))) null.value <- c(null.value,which(dat1[,i]<0))
null.value <- c(null.value,which(apply(dat1[,2:ncol(dat1)],1,sum,na.rm=T)==0))

if(length(null.value)!=0) {LogEvent.removed <- dat1[sort(unique(null.value)),] ; dat1 <- dat1[-sort(unique(null.value)),]}
Store(LogEvent.removed)

names(dat1)[-1]=unlist(lapply(strsplit(names(dat1[,-1]),"_"),function(x) x[[3]]))

#removing miscellaneous species

dat1 <- dat1[,!names(dat1)=="MZZ"]
save(dat1, file="dat1_2007.Rdata")
#load("dat1_2007.Rdata")

#-----------------------------
# II. EXPLORING THE VARIOUS METHODS FOR IDENTIFYING MAIN SPECIES AND KEEPING THEM IN THE DATA SET (STEP 1)
#-----------------------------
#EXPLORATION
explo=selectMainSpecies(dat1,analysisName,RunHAC=TRUE,DiagFlag=FALSE)

# Step 1 : selection of main species
Step1=extractTableMainSpecies(dat1,explo$NamesMainSpeciesHAC,paramTotal=95,paramLogevent=100)

save(explo,Step1,file="Explo_Step1.Rdata")

load("Explo_Step1.Rdata")

#-----------------------------
# III. STEP 2 - PCA - NO-PCA
#-----------------------------

# Step 2 : PCA

#for (option_step2 in c("PCA_70","PCA_SC","NO_PCA")) {
option_step2="PCA_70"

setwd(paste(path,analysisName,sep="/"))
if (!file.exists(option_step2)) dir.create(option_step2)
setwd(paste(path,analysisName,option_step2,sep="/"))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)      


if (option_step2=="PCA_70") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="70percents") else    # criterion="70percents"
if (option_step2=="PCA_SC") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="screetest") else    # criterion="screetest"
if (option_step2=="NO_PCA") Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="nopca",criterion=NULL)    

save(Step2,file="Step2.Rdata")

load("Step2.Rdata")

#-----------------------------
# IV. STEP 3 - CLUSTERING METHOD : HAC, CLARA OR KMEANS
#-----------------------------
#for (option_step3 in c("CLARA","KMEANS")) {
option_step3="CLARA"

setwd(paste(path,analysisName,option_step2,sep="/"))
if (!file.exists(option_step3)) dir.create(option_step3)
setwd(paste(path,analysisName,option_step2,option_step3,sep="/"))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)      


if (option_step3=="HAC")    Step3=getMetierClusters(Step1,Step2,analysisName,methMetier="hac",param3="euclidean",param4="ward") else    
if (option_step3=="CLARA")  Step3=getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="clara",param3="euclidean",param4=NULL) else    
if (option_step3=="KMEANS") Step3=getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="kmeans",param3=NULL,param4=NULL)    

save(Step3,file="Step3.Rdata")
#load("Step3.Rdata")
#} # end of step 3
#} # end of step 2





#-----------------------------
# V. STEP 4 - MERGING BACK TO EFLALO
#-----------------------------

#choosing the final option
setwd(paste(path,analysisName,sep=""))

load("PCA_70/CLARA/Step3.Rdata")

 if(!nrow(dat1)==nrow(Step3$LE_ID_clust)) print("--error : number of lines in step 3 not equal to input eflalo, please check!!--")

dat1 <- cbind(dat1,CLUSTER=Step3$LE_ID_clust[,"clust"])
#now reload the full data set

eflalo_ori[-sort(unique(null.value)),"CLUSTER"] <- Step3$LE_ID_clust[,"clust"]




#-----------------------------------------------------
# VI. STEP 5 - COMPARISON WITH ORDINATION METHODS
#-----------------------------------------------------

# load previous R objects (Step1,Step2,Step3)
setwd(paste(path,analysisName,sep="/"))
load("dat1_2007.Rdata")
load("Explo_Step1.Rdata")
option_step2="PCA_70"
setwd(paste(path,analysisName,option_step2,sep="/"))
load("Step2.Rdata")
option_step3="CLARA"
setwd(paste(path,analysisName,option_step2,option_step3,sep="/"))
load("Step3.Rdata")
compOrdin="CompOrdin"
if (!file.exists(compOrdin)) dir.create(compOrdin)
setwd(paste(path,analysisName,option_step2,option_step3,compOrdin,sep="/"))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)  
compMetiers=compareMultivariateClassificationToOrdination(dat=dat1,Step2=Step2,clusters=Step3$clusters$clustering)
save(compMetiers,file="compMetiers.Rdata")






#-----------------------------
# VII. STEP 6 - Predicting Metier for current year using clustering performed previous year
#-----------------------------


# load previous R objects (Step1,Step2,Step3)
setwd(path)
option_step2="PCA_70"
option_step3="CLARA"
load(paste(analysisName,"/",option_step2,"/",option_step3,"/Step3.Rdata",sep=""))
clust2007=Step3$clusters$clustering
rm(Step3)
gc()
#load(paste(analysisName,"/",option_step2,"/Step2.Rdata",sep=""))
load(paste(analysisName,"/Explo_Step1.Rdata",sep=""))
rm(explo)
gc()

# load your new dataset (called dat1 here)
 #! KEEPING ONLY LE_ID AND THE OUTPUT YOU WANT TO GET  (KG/EURO)
# load your own dataset (called dat1 here)
new_country <- "All"
new_year <- 2008
new_AreaCodename <- "3a4"
new_Gear <- c("OTB")
#load("All_eflalo_2008OTB3a4.Rdata")
load("dat1_2008.Rdata")
datPred=dat1
rm(dat1)
gc()


datPred <- dat1[,c("LE_ID",grep("EURO",names(dat1),value=T))]
rm(dat1)
gc()
datPred[is.na(datPred)]=0

#removing negative and null values
null.value <- vector()
for (i in grep("EURO",names(datPred))) null.value <- c(null.value,which(datPred[,i]<0))
null.value <- c(null.value,which(apply(datPred[,2:ncol(datPred)],1,sum,na.rm=T)==0))

if(length(null.value)!=0) {LogEvent.removed <- datPred[sort(unique(null.value)),] ; datPred <- datPred[-sort(unique(null.value)),]}
#Store(LogEvent.removed)

names(datPred)[-1]=unlist(lapply(strsplit(names(datPred[,-1]),"_"),function(x) x[[3]]))

#removing miscellaneous species

datPred <- datPred[,!names(datPred)=="MZZ"]

#create a new folder for step6
option_step5 = "Predict_2008"
setwd(paste(path,analysisName,option_step2,option_step3,sep="/"))
if (!file.exists(option_step5)) dir.create(option_step5)
setwd(paste(path,analysisName,option_step2,option_step3,option_step5,sep="/"))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)

dim(datPred)
class(datPred)
colnames(datPred)
le_id_datPred <- datPred[,"LE_ID"]
nbSpeciesDatPred <- ncol(datPred)-1
datPred <- datPred[,-1]
nameDatPredSpecies <- colnames(datPred)

datPred=as.matrix(datPred)
#datPred <- matrix(datPred,ncol=nbSpeciesDatPred,nrow=length(le_id_datPred))
rownames(datPred) <- le_id_datPred
#colnames(datPred) <- nameDatPredSpecies
Donnees2008Cluster = predictMetier(learningData=Step1,clustersAffectation=clust2007,newData=datPred)
save(Donnees2008Cluster,file="Donnees2008Cluster.Rdata")







sam=sample(1:nbLog,size=15000,replace=F) 
outofsam=setdiff(1:nbLog,sam)
learningData=dat[sam,]
learningData=learningData[,-1]
newData=as.data.frame(dat[outofsam,])
newData=newData[,-1]
clustersAffectation=Step3$clusters$clustering[sam]
learningData=as.data.frame(cbind(learningData,as.factor(clustersAffectation)))
colnames(learningData)=c(colnames(learningData[1:nbSpeciesLearningData]),"clust")

learning=nm(clust~.,data=learningData)      # celui là marche !
gam=learning$regularization["gamma"]
lam=learning$regularization["lambda"]
learning=rda(clust~.,data=learningData,gamma=gam,lambda=lam)
result=predict(learning,newData) 
clustPred=result$class
clustReal=Step3$clusters$clustering[outofsam]
t=table(clustPred,clustReal)
t=t[order(as.integer(rownames(t))),]
pourcWellClass=sum(diag(t))/length(outofsam)*100
pourcBadClass=100-pourcWellClass


learning=fda(clust~.,data=learningData)
result=predict(learning,newData) 
clustPred=result$class
clustReal=Step3$clusters$clustering[outofsam]
table(clustPred,clustReal)

model=multinom(clust~.,data=learningData)
stepc=greedy.wilks(clust~.,data=learningData,niveau=0.1)