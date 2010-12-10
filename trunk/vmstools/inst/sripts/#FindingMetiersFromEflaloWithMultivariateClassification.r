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

path <- "C:/CLARA/VMSTools_sandbox/"
setwd(path) # you must choose the path of your working directory

source("programs/FunctionsForClassif.r")
source("programs/Classif.r")
source("programs/ExploSpeciesSelection.r")
memory.limit(size=4000)



#-----------------------------
# I. GETTING THE DATA IN AND CLEANING FOR MISSING AND NEGATIVE VALUES ETC
#-----------------------------
 
country <- "SCO"

year <- 2007
AreaCodename <- "3a4"
Gear <- c("OTB")

analysisName=paste(country,"_",Gear,year,"_",AreaCodename,sep="")

# load your own dataset (called dat1 here)
load(paste("data/",country,"_eflalo_",year,Gear,AreaCodename,".Rdata",sep=""))


# combined data
dat <- read.table("data/FRA-DK-NLD_EFLALO_2007_OTB_EURO_NS-2.txt",dec=",",sep=";",header=TRUE)
dat$LE_ID <- paste(dat$COUNTRY,dat$LE_ID,sep=".")
dat <- dat[,c("LE_ID",grep("EURO",names(dat),value=T))]


sco <- read.csv("data/OTB_2007_SCO.csv")
sco$LE_ID <- paste("SCO",sco$LE_ID,sep=".")
sco <- sco[,c("LE_ID",grep("EURO",names(sco),value=T))]

combined_columns <- unique(c(names(dat),names(sco)))

dat[,combined_columns[is.na(match(combined_columns,names(dat)))]] <- NA
sco[,combined_columns[is.na(match(combined_columns,names(sco)))]] <- NA

dat1 <- rbind(dat,sco)

rm(dat,sco)
gc()


# creating the directory of the analysis
if (!file.exists(analysisName)) dir.create(analysisName)
setwd(paste(path,analysisName,sep=""))
#delete old R cache
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)      


eflalo_ori <- dat1 # keeping this in cached memory for making the final merging at the end
Store(eflalo_ori)

#! KEEPING ONLY LE_ID AND THE OUTPUT YOU WANT TO GET  (KG/EURO)
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


#-----------------------------
# II. EXPLORING THE VARIOUS METHODS FOR IDENTIFYING MAIN SPECIES AND KEEPING THEM IN THE DATA SET (STEP 1)
#-----------------------------

#EXPLORATION
explo=ExploSpeciesSelection(dat1,analysisName,RunHAC=TRUE,DiagFlag=FALSE)

# Step 1 : selection of main species
Step1=classif_step1(dat1,explo$NamesMainSpeciesHAC,paramTotal=95,paramLogevent=100)

save(explo,Step1,file="Explo_Step1.Rdata")


#-----------------------------
# III. STEP 2 - PCA - NO-PCA
#-----------------------------

# Step 2 : PCA

for (option_step2 in c("PCA_70","PCA_SC","NO_PCA")) {

setwd(paste(path,analysisName,sep=""))
if (!file.exists(option_step2)) dir.create(option_step2)
setwd(paste(path,analysisName,"/",option_step2,sep=""))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)      


if (option_step2=="PCA_70") Step2=classif_step2(Step1,analysisName,pcaYesNo="pca",criterion="70percents") else    # criterion="70percents"
if (option_step2=="PCA_SC") Step2=classif_step2(Step1,analysisName,pcaYesNo="pca",criterion="screetest") else    # criterion="screetest"
if (option_step2=="NO_PCA") Step2=classif_step2(Step1,analysisName,pcaYesNo="nopca",criterion=NULL)    

save(Step2,file="Step2.Rdata")


#-----------------------------
# IV. STEP 3 - CLUSTERING METHOD : HAC, CLARA OR KMEANS
#-----------------------------

for (option_step3 in c("CLARA","KMEANS")) {

setwd(paste(path,analysisName,"/",option_step2,sep=""))
if (!file.exists(option_step3)) dir.create(option_step3)
setwd(paste(path,analysisName,"/",option_step2,"/",option_step3,sep=""))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)      


if (option_step3=="HAC")    Step3=classif_step3(Step1,Step2,analysisName,methMetier="hac",param3="euclidean",param4="ward") else    
if (option_step3=="CLARA")  Step3=classif_step3(Step1,Step2,analysisName=analysisName,methMetier="clara",param3="euclidean",param4=NULL) else    
if (option_step3=="KMEANS") Step3=classif_step3(Step1,Step2,analysisName=analysisName,methMetier="kmeans",param3=NULL,param4=NULL)    

save(Step3,file="Step3.Rdata")
} # end of step 3
} # end of step 2




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




#-----------------------------
# VI. STEP 5 - Predicting Metier for current year using clustering performed previous year
#-----------------------------
 
new_country <- "SCO"
new_year <- 2008
new_AreaCodename <- "3a4"
new_Gear <- c("OTB")
# load previous R objects (Step1,Step2,Step3)
setwd(path)
load(paste(analysisName,"/PCA_70/CLARA/Step3.Rdata",sep=""))
load(paste(analysisName,"/PCA_70/Step2.Rdata",sep=""))
load(paste(analysisName,"/Explo_Step1.Rdata",sep=""))

# load your new dataset (called dat1 here)
#! KEEPING ONLY LE_ID AND THE OUTPUT YOU WANT TO GET  (KG/EURO)
load(paste("data/",new_country,"_eflalo_",new_year,new_Gear,new_AreaCodename,".Rdata",sep=""))
datPred <- dat1[,c("LE_ID",grep("EURO",names(dat1),value=T))]
datPred[is.na(datPred)]=0

#removing negative and null values
null.value <- vector()
for (i in grep("EURO",names(datPred))) null.value <- c(null.value,which(datPred[,i]<0))
null.value <- c(null.value,which(apply(datPred[,2:ncol(datPred)],1,sum,na.rm=T)==0))

if(length(null.value)!=0) {LogEvent.removed <- datPred[sort(unique(null.value)),] ; datPred <- datPred[-sort(unique(null.value)),]}
Store(LogEvent.removed)

names(datPred)[-1]=unlist(lapply(strsplit(names(datPred[,-1]),"_"),function(x) x[[3]]))

#removing miscellaneous species

datPred <- datPred[,!names(datPred)=="MZZ"]

#create a new folder for step5
option_step2="PCA_70"
option_step3="CLARA"
option_step5 = "Predict_2008"
setwd(paste(path,analysisName,"/",option_step2,"/",option_step3,sep=""))
if (!file.exists(option_step5)) dir.create(option_step5)
setwd(paste(path,analysisName,"/",option_step2,"/",option_step3,"/",option_step5,sep=""))
if (file.exists(".R_Cache")) unlink(".R_Cache",recursive=TRUE)


Donnees2008Cluster = metierPredict(learningData=Step1,clustersAffectation=Step3$clusters$clustering,newData=datPred)
