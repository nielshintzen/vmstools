# Practical 2i. Metier Identification in vmstools #
## Introduction ##
The European Common Fisheries Policy (CFP) calls for the implementation of an ecosystem-based approach to fisheries management (EBFM), to achieve good environmental status as required by the Marine Strategy Framework Directive. This has led in the last years to the enforcement of an area-based management scheme and includes measures to control and limit the fishing impact on the ecosystem. As part of the process, the CFP recognizes the importance of accounting for heterogeneity in fishing practices and there is a desire to move toward a fleet-based management. To make this approach operational, detailed information is needed on the various fishing activities of the fleets, as well as on their spatial and temporal distribution.  The overall aim is to pigeonhole fishing activity into fewer and more easy-to-manage categories. However, no single and fully unified approach has emerged, as there are so many criteria and scales that can be used. The European Community’s Data Collection Framework (DCF; EC, 2008), has been an important milestone and suggests a common framework for all European Member states, defining the needs and requirement for fisheries data collection. The DCF has adopted the following definitions which will be used throughout the subsequent analyses:

A fleet segment is a group of vessels with the same length class (LOA, Length Over All) and predominant fishing gear during the year. Vessels may have different fishing activities during the reference period, but might be classified in only one fleet segment.

A métier is a group of fishing operations targeting a similar (assemblage of) species, using similar gear, during the same period of the year and/or within the same area, and which are characterized by a similar exploitation pattern. [Note: ICES (2003) used to define the term fishery as well, as being the metier practiced by a given fleet segment, but this distinction has not been maintained].

The analysis of catch and effort per métier would allow more accurate estimates of fishing mortalities and so the DCF has thus defined metiers according to a hierarchical structure using six nested levels: level 1- Activity (fishing/non fishing);
level 2- Gear class (e.g., trawls, dredges, etc.);
level 3- Gear groups (e.g., bottom trawl, pelagic trawl);
level 4- EU level gear type (e.g., Bottom Otter trawl OTB, bottom pair trawl PTB);
level 5- Target assemblage (e.g., Demersal fish, crustaceans);
level 6- Selectivity feature of the gear (mesh size class and presence/absence of a selective device with its technical characteristics).

Here we define a seventh level which is not an official DCF definition but it represents the 'true' targeted species better than the rather crude DCF level 5.  To do this it is clearly necessary to establish common methods and tools for grouping logbook data into metiers. Logbook data provide information for each trip by each vessel on the amount of the main species caught and kept on board; the catch date; the specific ICES statistical rectangle; and the type of gear used.  The information on DCF levels 1 to 4, and to some extent level 6, is straightforward because it is directly available in the logbooks. However, fishermen do not declare which species they are actually targeting when fishing, and therefore the definition of Level 5, is expected to (but may not) reflect the fishing intention, and must be inferred from the available data, and in particular the catch data.

The simplest approach consists of selecting the fishing trips where a certain catch proportion of selected key species is exceeded. This approach is largely based on trial & error processes, and requires often a qualitative a priori knowledge on the fisheries. Another approach consists of conducting multivariate analyses on the species composition in catch data by trip or fishing operation (referred to as catch profiles), then grouping similar catch profiles into métiers. This grouping can be performed by direct visual inspection (Biseau and Gondeaux, 1988; Laurec et al., 1991) or statistically through clustering analysis.

Included in vmstools is an operational tool that: i) does sensitivity analyses of métier definition (at the level 7) to a variety of previously used classification methods; ii) predicts the métier of a new dataset; and iii) links the results obtained at level 7 to the DCF level 5.


## Preparation and cleaning ##

Make sure there is no rubbish lying around in R:

```

rm(list=ls())
```

Make sure libraries are attached:
```

library(vmstools)
library(maps)
library(mapdata)
```

Increase memory size to maximum:
```

memory.size(4000)
```


Set up directories where the code will leave the output:
```

codePath <- "C:/R/"
dataPath <- "C:/data/"
outPath  <- "C:/output/"
```

Set various parameters such as the maximum mesh size:

```

msThres  <- 150  #Maximum meshsize threshold in mm
lanThres <- 1.5  #Maximum difference in log10-transformed sorted weights
```

Make sure the eflalo data are clean:
```

data(eflalo)
idx <- which(eflalo$LE_MSZ >= msThres)
eflalo <- eflalo[-idx,]
```

Remove outlying catch records:
```

specs <- substr(colnames(eflalo)[grep("KG",colnames(eflalo))],7,9)
#- Calculate the landing boundaries
specBounds <- lapply(as.list(specs),function(x){ idx <- grep(x,colnames(eflalo))[grep("KG",colnames(eflalo)[grep(x,colnames(eflalo))])];
wgh <- sort(unique(eflalo[which(eflalo[,idx]>0),idx]));
difw <- diff(log10(wgh));
return(ifelse(any(difw > lanThres),wgh[rev(which(difw <= lanThres))],ifelse(length(wgh)==0,0,max(wgh,na.rm=T))))})

specBounds <- cbind(specs,unlist(specBounds));
specBounds[which(is.na(specBounds[,2])==T),2] <- "0"
#- Get the index of each of the species
idx <- unlist(lapply(as.list(specs),function(x){ idx   <- grep(x,colnames(eflalo))[grep("KG",colnames(eflalo)[grep(x,colnames(eflalo))])];
return(idx)}))
```

If landing is greater than threshold boundary then it gets to be NA:

```

for(iSpec in idx) eflalo[which(eflalo[,iSpec] > an(specBounds[(iSpec-idx[1]+1),2])),iSpec] <- 0
for(i in kgeur(colnames(eflalo))) eflalo[which(is.na(eflalo[,i]) == T),i] <- 0
```

Remove non-unique trip numbers:
```

eflalo <- eflalo[!duplicated(eflalo$LE_ID),]
```

Remove records with arrival date before departure date:
```

eflalop <- eflalo
eflalop$FT_DDATIM <- as.POSIXct(paste(eflalo$FT_DDAT,  eflalo$FT_DTIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S")
eflalop$FT_LDATIM <- as.POSIXct(paste(eflalo$FT_LDAT,  eflalo$FT_LTIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S")
idx  <- which(eflalop$FT_LDATIM >= eflalop$FT_DDATIM)
eflalo <- eflalo[idx,]
```


## Perform the métier analyses ##

Start by attaching inbuilt datasets to R search path:
```

data(correspLevel7to5)
data(correspMixedMetier)
```
Note: Have a look at these data to see if you can work out what they might do.

```

eflaloOrig <- eflalo # Good idea to 'save' the original data in case of  some mistake.
```

Give the analysis a name and create a directory for the output:

```

analysisName="Atlantis"
dir.create(outPath);
setwd(outPath)
dir.create(analysisName)
setwd(paste(outPath,analysisName,sep="/"))
```

Keep only L\_ID and the value columns for economy & tidiness:

```

eflalo <- eflalo[,c("LE_ID",grep("EURO",names(eflalo),value=T))]
eflalo[is.na(eflalo)]=0
```

Removing any negative and null values:
```

null.value <- vector()
for (i in grep("EURO",names(eflalo))) null.value <- c(null.value,which(eflalo[,i]<0))
null.value <- c(null.value,which(apply(eflalo[,2:ncol(eflalo)],1,sum,na.rm=T)==0))

if(length(null.value)!=0) {
LogEvent.removed  <- eflalo[sort(unique(null.value)),]
eflalo <- eflalo[-sort(unique(null.value)),]
Store(LogEvent.removed)
}
names(eflalo)[-1]=unlist(lapply(strsplit(names(eflalo[,-1]),"_"),function(x) x[[3]]))
```

Removing 'miscellaneous' species:
```

eflalo  <- eflalo[,!names(eflalo)=="MZZ"]
```

Optional save of the cleaned eflalo data to disk if you don't want to keep doing the above!
```

save(eflalo,file=paste(outPath,"cleanEflalo.RData",sep=""))
```


STEP 1 - Identifying main species

Exploration:
```

explo <-selectMainSpecies(eflalo,analysisName,RunHAC=TRUE,DiagFlag=FALSE)
```
STEP 1- Selection of main species
```

Step1<- extractTableMainSpecies(eflalo,explo$namesMainSpeciesHAC,paramTotal=95,paramLogevent=100)
```

Optional save of the Explo Step 1 output
```

save(explo,Step1,file="Explo_Step1.Rdata")
```

STEP 2 - Principal Components Analysis (PCA) or Not PCA

```

option_step2 <- "PCA_70"
if (!file.exists(option_step2)) dir.create(option_step2)
setwd(paste(outPath,analysisName,option_step2,sep="/"))

if (option_step2=="PCA_70") Step2 <- getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="70percents") else    # criterion="70percents"
if (option_step2=="PCA_SC") Step2 <- getTableAfterPCA(Step1,analysisName,pcaYesNo="pca",criterion="screetest") else    # criterion="screetest"
if (option_step2=="NO_PCA") Step2 <- getTableAfterPCA(Step1,analysisName,pcaYesNo="nopca",criterion=NULL)
```
Optional save of step 2 output:
```

save(Step2,file="Step2.Rdata")
```

STEP 3 - Clustering methods: HAC, CLARA or KMEANS (step 3)

Chose which clustering method to play with:
```

option_step3="CLARA"
```
Create directory to store output:
```

if (!file.exists(option_step3)) dir.create(option_step3)
setwd(paste(outPath,analysisName,option_step2,option_step3,sep="/"))
```

Do the cluster analyses:
```

if (option_step3=="HAC")    Step3 <- getMetierClusters(Step1,Step2,analysisName,methMetier="hac",param1="euclidean",param2="ward") else
if (option_step3=="CLARA")  Step3 <- getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="clara",param1="euclidean",param2=NULL) else
if (option_step3=="KMEANS") Step3 <- getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="kmeans",param1=NULL,param2=NULL)
```
Optional save of the step 3 data
```

save(Step3,file="Step3.Rdata")
```

Restore the data
```

eflalo <- eflaloOrig
```

## Exercises ##

Explore and interpret the output of the 3 steps.