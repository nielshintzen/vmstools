## Load in priority
rm(list=ls(all=TRUE))
gc(reset=TRUE)
#setwd("C:/CLARA/VMSTools_sandbox/")
setwd("Z:/VMSTools_sandbox/combinedFR_DK")    
source("../FunctionsForClassif.r")
source("../Classif.r")
#source("../ExploSpeciesSelection_clu.r")
memory.limit(size=2000)
#
#

## Load your data
year <- 2008
AreaCodename <- "3a4"
Gear <- c("OTB")
load(paste("eflalo_",year,Gear,AreaCodename,".Rdata",sep=""))
#which(dat1[,-1]<0)
#dat1 <- dat1[!is.na(dat1[,-1]),]
efl_fr <- read.table("EFLALO_2008_OTB_EURO_NS.txt",sep=";",header=TRUE)
efl_fr$VE_REF <- as.factor(efl_fr$VE_REF)
efl_fr$FT_REF <- as.factor(efl_fr$FT_REF)
efl_fr$LE_ID <- as.factor(efl_fr$LE_ID)


dat1 <- dat1[,c("VE_REF","FT_REF","LE_ID",grep("EURO",names(dat1),value=T))]

dat1[,names(efl_fr)[is.na(match(names(efl_fr),names(dat1)))]] <- NA
efl_fr[,names(dat1)[is.na(match(names(dat1),names(efl_fr)))]] <- NA

efl_comb <- rbind(dat1,efl_fr)
#efl_comb <- efl_comb[!is.na(efl_comb[,-1]),]

#removing negative values
null.value <- vector()
for (i in grep("EURO",names(efl_comb))) null.value <- c(null.value,which(efl_comb[,i]<0))
efl_comb <- efl_comb[-unique(null.value),]


rm(dat1,efl_fr)

gc(reset=TRUE)

#

t_init <- Sys.time()

#first exploring the method for species selection
#explo_species <- ExploSpeciesSelection(dat=dat1,analysisName="OTB2009",Val="EURO")

#explo_species <- ExploSpeciesSelection(dat=efl_comb,Val="EURO")
#with too many data the function doesn't work so we need to source a script instead
source("Z:/VMSTools_sandbox/ExploSpeciesSelection_Script.r")

gc(reset=TRUE)

#then choosing the preferred options and running step 1
load("explo_species.Rdata")
step1 <- classif_step1(dat=efl_comb,NamesMainSpeciesHAC=explo_species$NamesMainSpeciesHAC,Val="EURO",
                          paramTotal=95,paramLogevent=100)
                          
                          
save(step1, file="step1_2008OTB3a4.Rdata")
rm(explo_species)
Store(efl_comb)
gc(reset=TRUE)



#then exploring various options for step2

#load("step1_2008OTB3a4.Rdata")
analysisName = paste("pca70_hac_",year,sep="")
assign(analysisName,classif_step2(step1,analysisName,pcaYesNo="pca",
       criterion="70percents",methMetier="hac",param3="euclidian",param4="ward"))
save(list=analysisName,file=paste(analysisName,".Rdata",sep=""))
rm(list=analysisName)



analysisName = paste("pca70_kmeans_",year,sep="")
assign(analysisName,classif_step2(step1,analysisName,
       pcaYesNo="pca",criterion="70percents",methMetier="kmeans",param3=NULL,param4=NULL))
save(list=analysisName,file=paste(analysisName,".Rdata",sep=""))
rm(list=analysisName)
       
#
analysisName = paste("pca70_clara_",year,sep="")
assign(analysisName,classif_step2(step1,analysisName,pcaYesNo="pca",
       criterion="70percents",methMetier="clara",param3="euclidian",param4=NULL))
save(list=analysisName,file=paste(analysisName,".Rdata",sep=""))
rm(list=analysisName)
       
       

analysisName = paste("pcaScree_kmeans_",year,sep="")
assign(analysisName,classif_step2(step1,analysisName,pcaYesNo="pca",
       criterion="screetest",methMetier="kmeans",param3=NULL,param4=NULL))
save(list=analysisName,file=paste(analysisName,".Rdata",sep=""))
rm(list=analysisName)
       

analysisName = paste("pcaScree_clara_",year,sep="")
assign(analysisName, classif_step2(step1,analysisName,pcaYesNo="pca",
       criterion="screetest",methMetier="clara",param3="euclidian",param4=NULL))
save(list=analysisName,file=paste(analysisName,".Rdata",sep=""))
rm(list=analysisName)
       

analysisName = paste("pcaScree_hac_",year,sep="")
assign(analysisName,classif_step2(step1,analysisName,pcaYesNo="pca",
       criterion="screetest",methMetier="hac",param3="euclidian",param4="ward"))
save(list=analysisName,file=paste(analysisName,".Rdata",sep=""))
rm(list=analysisName)
       


print(Sys.time()-t_init)


