###################################################################################
##     FINDING METIERS FROM EFLALO USING MULTIVARIATE CLUSTERING METHODS         ##
##                                                                               ##
##              LOT 2 - VMS LOGBOOKS (WP2)                                       ##
##                                                                               ##
##                                                                               ##
## Authors : Nicolas Deporte (IFREMER / OCEANIC DEVELOPPEMENT)                   ## 
##           Stéphanie Mahévas, Sébastien Demanèche (IFREMER)                    ##
##           Clara Ulrich, Francois Bastardie (DTU Aqua)                         ##
##                                                                               ##
## Last update : February 2011   ,                                               ##
##                                                                               ##
## Runs with R 2.11.1                                                            ##
##                                                                               ##
###################################################################################



rm(list=ls(all=TRUE))
gc(reset=TRUE)

path <- "C:/Nicolas/Scripts/R/Analyses"
setwd(path) # you must choose the path of your working directory

source("multivariateClassificationUtils.r")
source("extractTableMainSpecies.r")
source("getTableAfterPCA.r")
source("getMetierClusters.r")
source("selectMainSpecies.r")
source("compareToOrdination.r")
source("predictMetier.r")
source("getEflaloMetierLevel7.r")
memory.limit(size=4000)

# Load the table linking 3A-CODE (FAO CODE of species) to the species assemblage (level 5).
load("correspLevel7to5.rda")
# Load the table linking mixed metiers (composed by 2 simple metiers) to their official code of mixed metiers level 5 (FAO - 3 characters).
load("correspMixedMetier.rda")

path <- "C:/Nicolas/Scripts/R/Analyses/Donnees_completes"
setwd(path)

#-----------------------------
# I. GETTING THE DATA IN AND CLEANING FOR MISSING AND NEGATIVE VALUES ETC
#-----------------------------

country <- "All"
year <- 2007
AreaCodename <- "3a4"
Gear <- c("OTB")

analysisName=paste(country,"_",Gear,year,"_",AreaCodename,sep="")

# Load your own dataset (called dat1 here)
load("All_eflalo_2007OTB3a4.Rdata")

eflalo_level7=getEflaloMetierLevel7(dat1, analysisName, path, critData="EURO", runHACinSpeciesSelection=TRUE, paramTotal=95, paramLogevent=100, critPca="PCA_70", algoClust="CLARA")

save(eflalo_level7,file="eflalo_level7.Rdata")
