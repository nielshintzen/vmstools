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


source("C:\\vmstools\\vmstools\\R\\multivariateClassificationUtils.r")
source("C:\\vmstools\\vmstools\\R\\selectMainSpecies.r")
source("C:\\vmstools\\vmstools\\R\\extractTableMainSpecies.r")
source("C:\\vmstools\\vmstools\\R\\getTableAfterPCA.r")
source("C:\\vmstools\\vmstools\\R\\getMetierClusters.r")
source("C:\\vmstools\\vmstools\\R\\compareMultivariateClassificationToOrdination.r")
source("C:\\vmstools\\vmstools\\R\\level7to5.r")
source("C:\\vmstools\\vmstools\\R\\level5.r")
source("C:\\vmstools\\vmstools\\R\\predictMetier.r")
memory.limit(size=4000)



#-----------------------------
# I. GETTING THE DATA IN AND CLEANING FOR MISSING AND NEGATIVE VALUES ETC
#-----------------------------
# 
#path <-"G:/DossiersASynchroniser/Stat-Methode-Mequapro/EtudePCP/R/RDec2010/"
path <- "C:/CLARA/VMSTools_sandbox/test Jan2011/"

country <- "All"
year <- 2007
AreaCodename <- "3a4"
Gear <- c("OTB")

analysisName=paste(country,"_",Gear,year,"_",AreaCodename,sep="")

setwd(file.path(path,analysisName) # you must choose the path of your working directory


# load your own dataset (called dat1 here)
#load(paste("data/All_eflalo_2007OTB3a4.Rdata",sep=""))
load("C:\\CLARA\\VMSTools_sandbox\\data\\All_eflalo_2007OTB3a4.Rdata")


eflalo_level7 <- getEflaloMetierLevel7(dat1,analysisName,path,critData="EURO", runHACinSpeciesSelection=FALSE, critPca="PCA_70", algoClust="CLARA")

save(eflalo_level7)
