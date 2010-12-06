rm(list=ls(all=TRUE))
gc(reset=TRUE)
setwd("C:/Analyses") # you must choose the path of your working directory
source("FunctionsForClassif.r")
source("Classif.r")
source("ExploSpeciesSelection.r")
memory.limit(size=4000)


dat=read.table(file="EFLALO_2008_DRB_EURO_NA.txt", sep=";",header = TRUE, quote="\"", dec=".")
dat=dat[,3:ncol(dat)]
dat[is.na(dat)]=0
names(dat)[-1]=unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))

analysisName="DRB"

# Exploration of main species
explo=ExploSpeciesSelection(dat,analysisName,RunHAC=TRUE,DiagFlag=FALSE)


# Step 1 : selection of main species
classif1=classif_step1(dat,explo$NamesMainSpeciesHAC,paramTotal=95,paramLogevent=100)


# Step 2 : PCA
classif2=classif_step2(classif1,analysisName,pcaYesNo="pca",criterion="70percents")    # criterion="70percents"
# or
classif2=classif_step2(classif1,analysisName,pcaYesNo="pca",criterion="screetest")    # criterion="screetest"


# Step 3 : clustering, and determination of metiers
classif31=classif_step3(classif1,classif2,analysisName,methMetier="hac",param3="euclidean",param4="ward")  # methMetier="hac"
# or
classif32=classif_step3(classif1,classif2,analysisName,methMetier="kmeans",param3=NULL,param4=NULL) # methMetier="kmeans"
# or
sam=sample(1:nrow(dat),5000,replace=F)
classif33=classif_step3(classif1[sam,],classif2[sam,],analysisName,methMetier="pam",param3=NULL,param4=NULL) # methMetier="pam"
# or
classif34=classif_step3(classif1,classif2,analysisName,methMetier="clara",param3="euclidean",param4=NULL) # methMetier="clara"
