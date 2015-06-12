# Practical 2: Cleaning and processing the VMS and logbook data #

## Introduction ##
Typically VMS data contain lots of errors and many of them have a technical origin. Some are easy to detect, others a bit more tricky. There are problems with logbook data too, but most often these are caused by the fishermen filling them out incorrectly, rather than by the processing and storage of the data.

Here we will provide some guidance on how to clean VMS (tacsat) and logbook (eflalo) data. It needs to be said, that these guidelines are not strict rules, and that other rules are most likely more appropriate for your own dataset!

## Cleaning VMS data ##
Again, we start by loading / reading in the VMS data:

```

data(tacsat)
```

In VMS we can distinguish roughly two types of errors, 1- those errors that are without doubt incorrect and 2- more ambiguous errors about which a bit of thinking may be required.

Among the errors that are incorrect, no matter what, are: VMS positions that are not actually on planet Earth; records that are true duplicates; and records that, either have headings outside a compass range, or are situated on land. Let's see if we can get rid of such erroneous records first.

```

idx <- which(abs(tacsat$SI_LATI) > 90 | abs(tacsat$SI_LONG) > 180) #points not on the globe
idx             <- unique(c(idx,which(tacsat$SI_HE < 0 | tacsat$SI_HE >  360))) #adding points with heading outside compass range
length(idx)
```

Luckily, no problems as yet. On to points on land and duplicate records, which needs a bit more code.

```

tacsat$SI_DATIM <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M") #create one date-time stamp
uniqueTacsat    <- paste(tacsat$VE_REF,tacsat$SI_LATI,tacsat$SI_LONG,tacsat$SI_DATIM) #get records as a string to easily check for duplicates
print(nrow(tacsat))

tacsat          <- tacsat[!duplicated(uniqueTacsat),] #get rid of the duplicates
print(nrow(tacsat))
```

And again, no problems, no duplicates (the only time you will see this most likely!). On to the points on land.

```

data(europa) #a dataset already build into the package with the coastlines embedded
head(europa)

#Convert the europa dataset into a polygon, takes rather long so be patient!
pols    <- lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),
function(x){data.frame(SI_LONG=subset(europa,SID==x)$X,
SI_LATI=subset(europa,SID==x)$Y)}))
idx     <- pointOnLand(tacsat,pols);
table(idx)

pol     <- tacsat[which(idx==1),] #points on land
#Overwrite the original tacsat dataset with only those points not on land
tacsat  <- tacsat[which(idx==0),] #tacsat not on land
```

Let's have a look where these points on land (according to our algorithm) are really located!

```

library(maps);library(mapdata)
map("worldHires",res=0,fill=T,col="darkgreen",xlim=c(-4,10),ylim=c(48,62)); map.axes()
points(x=pol$SI_LONG,y=pol$SI_LATI,col="red",pch=19,cex=0.5)
```

As you can see, most of them are either on land, close to a harbour or on a river. Hence, it is important to use a detailed map to do this assessment! You can specify your own map, you don't need to use the map supplied with VMStools.

This means we are half-way through cleaning the VMS dataset. Now we must investigate those types of 'errors' which are ambiguous and need to be the subject of further in-depth investigation and discussion. These are things like: VMS records associated with very high speeds; points in a harbour; and 'pseudo-duplicates', ie. points that are not exactly duplicates, but may be considered as such.

```

hist(tacsat$SI_SP,breaks=100) #are there any outlying results?
spThres <- 20 #I decide that 20knots is the maximum speed allowed in the dataset
idx <- which(tacsat$SI_SP > spThres)
tacsat <- tacsat[-idx,]
```

Very few points are associated with speeds higher than 20 knots. Let's move to the locations recorded in harbours. Obviously, these points are not true errors. Vessels are in harbour for large amounts of time. However, when analysing VMS data we are hardly ever interested in these points, and for that reason, we usually want to omit them. There are a large number of ports in Europe, checking all of them is not easy, but we will give it a try! Again, to make things simpler, a large dataset of harbour/port locations is already available as part of the VMStools package. Note, however, that this list is NOT complete, and that you might need to add additional harbours/ports (for your country) if you want to accurately perform this part. In particular the port dataset is only for the EU, and we were unable to get data for Portugal at the time.

```

data(euharbours)
idx <- pointInHarbour(tacsat$SI_LONG,tacsat$SI_LATI,harbours,saveHarbourList=F)
pih <- tacsat[which(idx==1),]
table(idx)
#Overwrite the tacsat file with only those records which are not in harbour
tacsat <- tacsat[which(idx==0),]
```

All in all, 5403 points are located in harbour (5% of the total!). If you've used the raw dataset, instead of the cleaned one, you will end up with 14787 rows instead. Do you know why?. Let's see where they really are located and if our analyses make sense.

```

library(maps);library(mapdata)
map("worldHires",res=0,fill=T,col="darkgreen",xlim=c(-4,10),ylim=c(48,62)); map.axes()
points(x=pih$SI_LONG,y=pih$SI_LATI,col="red",pch=19,cex=0.5)
```

### Exercise 1 ###
  1. Note the dot-at-sea in the corner with the Netherlands, Germany and  Denmark. Have a look in Google Maps to see which harbour is located there
  1. How many harbours are stored in the built-in harbour dataset?
  1. What does 'range' mean in the harbour dataset?

What remains in the task of cleaning the VMS dataset are the 'pseudo duplicates'. It is possible that the VMS transmitter on-board vessels can malfunction. Data might be transmitted at irregular intervals, e.g. 1, 2 or 4 hours. Also sometimes we see another position within a few seconds or few minutes of the first point. If you consider these to be 'true' errors, you might decide to remove them.
```

tacsat <- sortTacsat(tacsat) #sort the data by vessel and time
tacsatp<- intervalTacsat(tacsat,level="vessel",fill.na=T)

hist(tacsatp$INTV,breaks=1000000,xlim=c(0,150)) #apparently there is a large number of records with a very small interval rate (close to zero)

#I believe that interval rates lower than 5 minutes need to be removed
intThres <- 5
tacsat          <- tacsatp[which(tacsatp$INTV > intThres),]
```

Note that this cleaning describes the absolute basic parts you need to screen for.

### Exercise 2 ###
  1. Think of two other aspects you might need to screen for.
  1. Check your data for these aspects, and create appropriate graphs to show your result

## Cleaning logbook data ##
In cleaning the logbook data, similar issues as with cleaning the VMS data arise. Some errors are true mistakes and need to be fixed / removed, others depend on interpretation by the scientist. However, the origin of logbook errors are most often human mistakes rather than technical malfunctions as we have seen with VMS data. For this reason it is advisable to have a fundamental knowledge of the fisheries data and how these sources are collected and stored.

The basic cleaning of logbook data consists of removing records where arrival date occurs before departure date, duplicate records and easy to detect mistakes in the catch values recorded. Other aspects such as mesh sizes and vessel lengths or horse power being outside a 'believable' range can be dealt with as well.

```

data(eflalo)
#In principle, the LE_ID identifier should be unique, this allows for an easy solution
#Checking for duplicates
eflalo <- eflalo[!duplicated(eflalo$LE_ID),]

#If things are a bit more difficult
ID <- paste(eflalo$VE_REF,eflalo$VE_FLT,eflalo$VE_COU,eflalo$VE_LEN,eflalo$VE_KW,eflalo$VE_KW,eflalo$VE_TON,eflalo$FT_REF,eflalo$FT_DCOU,eflalo$FT_DHAR,eflalo$FT_DDAT,eflalo$FT_DTIME,eflalo$FT_LCOU,eflalo$FT_LHAR,eflalo$FT_LDAT,eflalo$FT_LTIME)
eflalo <- eflalo[!duplicated(ID),]
```

How about removing those records where arrival date occurs before departure date

```

eflalop           <- eflalo

#Create a date-time stamp for the departure date
eflalop$FT_DDATIM <- as.POSIXct(paste(eflalo$FT_DDAT,  eflalo$FT_DTIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M") #turn date and time into one date-time stamp

#Create a date-time stamp for the landing date
eflalop$FT_LDATIM <- as.POSIXct(paste(eflalo$FT_LDAT,  eflalo$FT_LTIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M") #turn date and time into one date-time stamp

#Now see where the landing date happens before departure
idx               <- which(eflalop$FT_LDATIM >= eflalop$FT_DDATIM)
print(nrow(eflalo))

#only keep the records we want
eflalo            <- eflalo[idx,]
print(nrow(eflalo))
```

Because we carefully screened the example Atlantis dataset before it was included in the package, you don't run into problems here.

Due to human mistakes, some catch records in the eflalo data end up being extremely large. Do we trust these values or not? Perhaps, an easy check is to see if these large catches are really big compared to other catches of the same species. If the second highest catch is more than 10 times smaller than the biggest catch you can start to question the results. Again, knowledge on the fishery and their landings is required!
As a rule of thumb, I normally use a factor 30 difference as a quick check.

```

#First get the species names in your eflalo dataset
specs  <- substr(colnames(eflalo[grep("KG",colnames(eflalo))]),7,9)
lanThres <- 1.5 #this is approx 30x difference on a log10 scale)

#Define per species what the maximum allowed catch is (larger than that value you expect it to be an error / outlier
specBounds      <- lapply(as.list(specs),function(x){
idx   <- grep(x,colnames(eflalo))[grep("KG",colnames(eflalo)[grep(x,colnames(eflalo))])];
wgh   <- sort(unique(eflalo[which(eflalo[,idx]>0),idx]));
difw  <- diff(log10(wgh));
return(ifelse(any(difw > lanThres),wgh[rev(which(difw <= lanThres)+1)],ifelse(length(wgh)==0,0,max(wgh,na.rm=T))))})

#Make a list of the species names and the cut-off points / error / outlier point
specBounds      <- cbind(specs,unlist(specBounds));

#Put these values to zero
specBounds[which(is.na(specBounds[,2])==T),2] <- "0"

#Get the index (column number) of each of the species
idx             <- unlist(lapply(as.list(specs),function(x){
idx   <- grep(x,colnames(eflalo))               [grep("KG",colnames(eflalo)[grep(x,colnames(eflalo))])];
return(idx)}))

#If landing > cut-off turn it into an 'NA'
for(iSpec in idx) eflalo[which(eflalo[,iSpec] > an(specBounds[(iSpec-idx[1]+1),2])),iSpec] <- NA

#Turn all other NA's in the eflalo dataset in KG and EURO columns to zero
for(i in kgeur(colnames(eflalo))) eflalo[which(is.na(eflalo[,i]) == T),i] <- 0
```

You have now finished cleaning the eflalo dataset, first save the results for later use, and time for an exercise.

```

save(tacsat,file="D:/output/tacsatClean.RData")
save(eflalo,file="D:/output/eflaloClean.RData")
```

### Exercise 3 ###
  1. Clean the eflalo dataset (remove these records) where meshsizes are larger than 150mm, how many do you come across?
  1. There are some very large vessels (> 100meters) in the eflalo dataset. What type of gear do they use and what do they primarily catch?