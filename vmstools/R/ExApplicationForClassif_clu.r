## Load in priority
rm(list=ls(all=TRUE))
gc(reset=TRUE)
#setwd("C:/CLARA/VMSTools_sandbox/")
setwd("Z:/VMSTools_sandbox/")    
source("FunctionsForClassif.r")
source("Classif.r")
source("ExploSpeciesSelection.r")
memory.limit(size=2000)
#
#

## Load your data
year <- 2009
AreaCodename <- "3a4"
Gear <- c("OTB")
load(paste("eflalo_",year,Gear,AreaCodename,".Rdata",sep=""))
#
## Compute the classification thanks to the function "classif"
## FIRSTLY, YOU HAVE TO READ THE BEGINNING OF THE SOURCE FILE "CLASSIF.R" TO KNOW HOW TO USE THE FUNCTION "CLASSIF"
## For example, if you have chosen HAC for the species selection, PCA and CLARA for the clustering of logevents, you have to write :

t_init <- Sys.time()

#first exploring the method for species selection
#explo_species <- ExploSpeciesSelection(dat=dat1,analysisName="OTB2009",Val="EURO")
explo_species <- ExploSpeciesSelection(dat=dat1,Val="EURO")


gc(reset=TRUE)

#then choosing the preferred options and running step 1

step1 <- classif_step1(dat=dat1,Val="EURO",methSpecies="logevent",param1=100,param2=NULL)
gc(reset=TRUE)


#then exploring various options for step2



pca70_kmeans <- classif_step2(step1,analysisName="PCA70_KM",pcaYesNo="pca",criterion="70percents",methMetier="kmeans",param3=NULL,param4=NULL)

pca70_pam <- classif_step2(step1,analysisName="PCA70_PAM",pcaYesNo="pca",criterion="70percents",methMetier="pam",param3=NULL,param4=NULL)

pca70_clara <- classif_step2(step1,analysisName="PCA70_CLARA",pcaYesNo="pca",criterion="70percents",methMetier="clara",param3="euclidian",param4=NULL)

pca70_hac <- classif_step2(step1,analysisName="PCA70_HAC",pcaYesNo="pca",criterion="70percents",methMetier="hac",param3="euclidian",param4="ward")


print(Sys.time()-t_init)


