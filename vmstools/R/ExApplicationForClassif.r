# Load in priority
rm(list=ls(all=TRUE))
gc()
setwd("C:/folder") # you must choose the path of your working directory
source("FunctionsForClassif.r")
source("Classif.r")
memory.limit(size=4000)


# Load your data
dateflalo=read.table(file="EflaloData.txt", sep=";",header = TRUE, quote="\"", dec=",")
dateflalo=dateflalo[,3:ncol(dateflalo)]
# Replacing NA's by 0
dateflalo[is.na(dateflalo)]=0
head(dateflalo)
dim(dateflalo)


# Compute the classification thanks to the function "classif"
# FIRSTLY, YOU HAVE TO READ THE BEGINNING OF THE SOURCE FILE "CLASSIF.R" TO KNOW HOW TO USE THE FUNCTION "CLASSIF"
# For example, if you have chosen HAC for the species selection, PCA and CLARA for the clustering of logevents, you have to write :
dat=dateflalo
analysisName=""
methSpecies="hac"
param1="euclidean"
param2="ward"
pcaYesNo="pca"
criterion="70percents"
methMetier="clara"
param3="euclidean"
param4=NULL
cahacpclara=classif(dat,analysisName,methSpecies,param1,param2,pcaYesNo,criterion,methMetier,param3,param4)


# Example
OTB2008euroNS=read.table(file="EFLALO_2008_OTB_EURO_NS.txt", sep=";",header = TRUE, quote="\"", dec=".")
OTB2008euroNS=OTB2008euroNS[,3:ncol(OTB2008euroNS)]
OTB2008euroNS[is.na(OTB2008euroNS)]=0
dat=OTB2008euroNS 
classif_OTB2008euroNS=classif(dat,"OTB2008euroNS","hac","euclidean","ward","pca","70percents","clara","euclidean",NULL)
summary(classif_OTB2008euroNS)