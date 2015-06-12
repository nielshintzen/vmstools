# Solutions of exercises #

## Practical 1 ##

### Exercise 1 ###
  1. See the tacsat and eflalo description (available under the download tab) for more info
  1. 'chr' stands for 'character', a string of letters. 'num' stands for 'numeric', which is a number with decimal places.
  1. Perhaps you've seen 'int' for 'integer' which is a full real number (no digits). You'll come across 'list', 'factor', 'matrix', 'array' or 'data.frame' a lot.

### Exercise 2 ###
  1. See the tacsat and eflalo description (available under the download tab) for more info
  1. The two new classes are: 'factor' and 'logi' for logical. Logical can only take two value, either TRUE or FALSE. Factor might look like a number or a character, but is neither of both. It is a sort of category classification, like in clinical tests where you might have a group A and B which have e.g. a treatment in common but both groups come from different cities.
  1. That would be 3050242 KG of plaice. You'll get that by using the following code: `sum(eflalo$LE_KG_PLE,na.rm=T)`

### Exercise 3 ###
  1. Run the `str(eflalo)` command to see if your formats are correct
  1. If columns out of your eflalo dataset are not yet in the right format, this function will convert those columns into the right format. Check before and after using the function with `str(eflalo)`.
  1. I hope so, but please consult one of the instructors if you get stuck

### Exercise 4 ###
1. This is a bit tricky, because the vessel names are 'characters' and therefore not numbers. R is not that picky howver but let's do it the right way here.

```

#Create a new column with numeric names for the vessel ID's
tacsat$numericName <- as.numeric(tacsat$VE_REF)
tacsatPartOne      <- subset(tacsat,VE_REF <  1000)
tacsatPartTwo      <- subset(tacsat,VE_REF >= 1000)

#Save to disk
#First set path where I want to save these files
path2Save          <- "D:/output/" #different on everyones computer
write.table(tacsatPartOne, file = file.path(path2Save,"tacsatPartOne.csv"), quote = TRUE, sep = ",",
eol = "\n", na = "NA", dec = ".", row.names = F,col.names = TRUE)
write.table(tacsatPartTwo, file = file.path(path2Save,"tacsatPartTwo.csv"), quote = TRUE, sep = ",",
eol = "\n", na = "NA", dec = ".", row.names = F,col.names = TRUE)
```

2. We read both files in using the readTacsat function and thereafter combine them

```

tacsatPartOne      <- readTacsat(file.path(path2Save,"tacsatPartOne.csv"))
tacsatPartTwo      <- readTacsat(file.path(path2Save,"tacsatPartTwo.csv"))

#Combine the two datasets
totalTacsat        <- rbindTacsat(tacsatPartOne,tacsatPartTwo)
```

3. Again 97015 rows of data with 9 columns. That is similar to: `data(tacsat); nrow(tacsat)`

## Practical 2 ##

### Exercise 1 ###
  1. That is Helgoland harbour
  1. In total 3839. You'll get to this number by checking: `data(euharbours); nrow(harbours)`
  1. Please read the help file of 'euharbours' `?euharbours`  The range defines the radius, from the harbour midpoint, the harbour spans. I.e. every VMS point inside this 'range' is classified as being in harbour.

### Exercise 2 ###
  1. Have a look at the Nations of vessels present in your dataset. Often you can have access to VMS pings of foreign vessels but not have access to logbook data. You could also screen for artefacts in GPS positions. I've seen for example simple plots of VMS points with strange patterns popping up because a mistake was made in converting from minute to decimal notation.
  1. `table(tacsat$VE_COU)`; `plot(tacsat$SI_LONG,tacsat$SI_LATI,pch=".",type="p")`

### Exercise 3 ###
1. This is how it's done:

```

#Find tacsat records with mesh size > 150mm
idx <- which(eflaloLE_MSZ > 150)
#Remove them from your file
eflalo <- eflalo[-idx,]
#count how many there were
print(length(idx))
```

2. Let's make a selection of these vessels and then look at their characteristics in two ways
```

#eflaloLV = eflalo Large Vessels
eflaloLV <- subset(eflalo,VE_LEN > 100)

#What are their characteristics
summarizeEflalo(eflaloLV)

#Find the gears
table(eflaloLV$LE_GEAR)

#Find their catch
sort(colSums(eflaloLV[,kgeur(colnames(eflaloLV))],na.rm=T))
```

## Practical 3 ##

### Exercise 1 ###
1. The main check is to look at the FT\_REF column that has just been added to the tacsatp dataset. It may contain some zero's but should also contain number > 0.
```

head(tacsatp) # I notice some FT_REF > 0 already
length(which(tacsatp$FT_REF != 0)) / nrow(tacsatp) #Apparently, nearly 97% of all tacsat records is associated with an FT_REF > 0 value

# Cross-check with eflalo
subT <- subset(tacsatp,FT_REF == tacsatp$FT_REF[1])
subE <- subset(eflalo[,1:25],FT_REF == tacsatp$FT_REF[1])
#Check vessel name
unique(subT$VE_REF) == unique(subE$VE_REF)
#-Check times in tacsat and eflalo
range(subT$SI_DATIM); subE[,c("FT_DDAT","FT_DTIME","FT_LDAT","FT_LTIME")]
```

2. Check with this line of code: `table(tacsatp$LE_GEAR)`

### Exercise 2 ###
  1. There are a number of differences, so it's not a surprise there is a large difference in KG's. First checking the KG's in both datasets by: `sum(eflalo$LE_KG_PLE)` and `sum(tacsatEflalo$LE_KG_PLE)`
> The difference in KG is approx. 1368879. This because the eflalo dataset spans both the years 1800 and 1801. Also, be aware that if in the 'splitAmongPings' function, the 'conserve' option would be turned off,
> the KGs would be even lower in the tacsatEflalo dataset.

## Practical 4 ##

### Exercise 1 ###
  1. It means that some of the rectangles could not be converted into longitude - latitude positions.
  1. The trick here is that converting a GPS position to an ICES rectangle is based on a simple calculation. Even if GPS positions
> are located outside the ICES region, the simple calculation can still be applied (but results in odd rectangle names). The conversion
> back however, relies on a bit more than a simple calculation and hence throws you a warning for the odd rectangle names out there.
> Below you'll find where these points are located
```

idx <- which(is.na(pos$SI_LONG)==T)
require(maps)
require(mapdata)
map("worldHires",xlim=range(tacsat$SI_LONG[idx],na.rm=T),ylim=range(tacsat$SI_LATI[idx]),fill=T,col="darkgreen")
points(tacsat$SI_LONG[idx],tacsat$SI_LATI[idx],col=2,pch=19,cex=0.3)
```

### Exercise 2 ###
1. The code to do this is rather similar to the aggregation of plaice in the example above. This is how it's done:
```

#Add 0.5 degree notation to tacsatEflalo
tacsatEflalo$LE_SPAT <- CSquare(tacsatEflalo$SI_LONG,tacsatEflalo$SI_LATI,degrees=0.5)
#Aggregate the COD catches
aggregate(tacsatEflalo$LE_KG_COD,by=list(tacsatEflalo$LE_SPAT),FUN=sum,na.rm=T)
```

2. Have a look at Practical 9, where a template for a plotting routine is given

### Exercise 3 ###
1. We first make a summary of all gridcells present in the tacsat dataset, list the highest value and thereafter find the lon-lat position.
```

#Get a summary of all gridcells present in the tacsat dataset
sort(table(overlay(Grid,SpatialPoints(coordinates(tacsat[,c("SI_LONG","SI_LATI")])))))

#Apparently, gridcell 94 has 7647 records associated
coordinates(Grid)[94,] # 4-longitude by 54 latitude
```

## Practical 5 ##

### Exercise 1 ###
1. Here is the code to plot the second interpolation.
```

#- Set the xlim and ylim of the plot
xrange <- range(c(cHs[[3]][-1,1],SL[[3]][-1,1]),na.rm=T)
yrange <- range(c(cHs[[3]][-1,2],SL[[3]][-1,2]),na.rm=T)
plot(cHs[[3]][-1,1],cHs[[3]][-1,2],type="l",xlim=xrange,ylim=yrange,xlab="Longitude",ylab="Latitude",
asp=1/lonLatRatio(mean(xrange),mean(yrange)))
lines(SL[[3]][-1,1],SL[[3]][-1,2],col=2,lty=2)
points(VMShf$SI_LONG[cHs[[3]][1,1]],VMShf$SI_LATI[cHs[[3]][1,1]],pch=19,col="blue")
points(VMShf$SI_LONG[cHs[[3]][1,2]],VMShf$SI_LATI[cHs[[3]][1,2]],pch=19,col="blue")
```
The interpolation is rather different because the heading at both VMS positions is quite different from the heading you would expect when
a vessel would steam straight to the next VMS point.

2. `distanceInterpolation(cHs[3])` or `distance(cHs[[3]][2:100,1],cHs[[3]][2:100,2],cHs[[3]][3:101,1],cHs[[3]][3:101,2])` and `distanceInterpolation(SL[3])`
3. The speed can be calculated if we know the time it took to get from point x to point y and if we know the distance travelled:
```

#Calculate the distance of the interpolation
distcHs <- distanceInterpolation(cHs[3])
distSL  <- distanceInterpolation(SL[3])
#Calculate the time difference between the two connected VMS pings (but first make a date-time column in the VMShf dataset
VMShf$SI_DATIM <- as.POSIXct(paste(VMShf$SI_DATE,VMShf$SI_TIME),format="%d/%m/%Y %H:%M")
timecHs <- difftime(VMShf$SI_DATIM[cHs[[3]][1,2]],VMShf$SI_DATIM[cHs[[3]][1,1]],units="mins")
#Calculate associated speed
distcHs / (as.numeric(timecHs)/60) #7.8 km/h
distSL  / (as.numeric(timecHs)/60) #6.2 km/h
```

## Practical 6 ##

### Exercise 1 ###
We take the approach as outlined above and follow these steps to get to the three questions under exercise 1.

  1. clean the tacsat and eflalo data
  1. merge them together
  1. select bottom trawling gear only
  1. perform activity analyses and take fishing activity only
  1. define the width of the gear
  1. define the grid

```

library(vmstools)
data(tacsat); data(eflalo)


# 1) We skip cleaning the dataset here, but load our pre-saved results
load("D:/output/tacsatClean.RData")
load("D:/output/eflaloClean.RData")

# 2) Merging eflalo and tacsat
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]

# 3) Select bottom gear only
tacsatp <- subset(tacsatp,LE_GEAR == "TBB" & year(tacsatp$SI_DATIM) == 1800)
eflalo  <- subset(eflalo,LE_GEAR == "TBB")

# 4) Perform activity
storeScheme <- activityTacsatAnalyse(tacsatp, units = "year", analyse.by = "LE_GEAR",
identify="means")

res           <- activityTacsat(tacsatp,units="year",analyse.by="LE_GEAR",storeScheme,
plot=TRUE,level="all")
tacsatp$SI_STATE <- res
# Take fishing pings only
tacsatp     <- subset(tacsatp,SI_STATE == "f")

# 5) Define with of the gear
# Often, you have a gear width column in your dataset. In that case, the line below
# adds the withd of a gear to the tacsat file.
tacsatp$LE_WIDTH <- eflalo$LE_WIDTH[match(tacsatp$FT_REF,eflalo$FT_REF)]

# In our case however, the eflalo dataset does not contain an LE_WIDTH column, so we have to make something up
tacsatp$LE_WIDTH <- 2*12 #two beams of each 12m wide

# 6) Define the grid
grd <- createGrid(xrange=range(tacsatp$SI_LONG),
yrange=range(tacsatp$SI_LATI),
resx=0.1,resy=0.05,type="SpatialGridDataFrame") #Let's take a tenth of an ICES rectangle for a start

# 7) Turn your tacsat points into SpatialPoints and project them on the grid
coords <- SpatialPoints(coordinates(tacsatp[,c("SI_LONG","SI_LATI")]))
tacsatp$grID <- over(coords,as(grd,"SpatialGrid"))

#8) Calculate area trawled per grid cell via the aggregate function
kmTrawled <- aggregate(tacsatp$INTV/60 * tacsatp$SI_SP * 1.852 * tacsatp$LE_WIDTH/1000,
by=list(tacsatp$grID),FUN=sum,na.rm=T)
colnames(kmTrawled) <- c("grID","km2")
```
1) Ok, we're there, that was a lot of work for calculation 1. The answer is: `sum(kmTrawled$km2)`

2) Now on to number 2, which is a bit simpler given our preparation above.
```

#9) Put the area trawled inside the grid data.frame
grd@data[kmTrawled$grID,1] <- kmTrawled$km2

#10) Calculate the surface area of each grid cell
grd  <- surface(grd)

#11) Now we can calculate proportion of total surface trawled (by gridcell)
grd@data$proptrawled <- grd@data$data / grd@data$cellArea
hist(grd@data$proptrawled[grd@data$proptrawled>0],breaks=100,xlab="Proportion trawled")
```

3) Finally, the 3rd question is relatively easy to get out of the SpatialGridDataFrame above
```

idx <- which(grd@data$proptrawled > 1)
print(length(idx))
```
So apparently, the anwer is 18, not a lot (given that the full grid contained nearly 11000 gridcells)

## Practical 7 ##

### Exercise 1 ###

The only part where we can supply a solution is on importing a shapefile.
As an example you can use the code below, but you need to adapt it to your own settings.

```

library(PBSmapping)
path        <- "D:/output/"
myShapeFile <- readShapePoly(paste(path,"myShapeFile.shp",sep=""))
plot(myShapeFile)
```

## Practical 8 ##

### Exercise 1 ###
A lot of the exercise is already performed above and we only need to make a small modification.
However, below I show the entire code workflow to arrive at the correct answer.

```

#Only use the VMS data in this instance
data(tacsat)
tacsat <- sortTacsat(tacsat)

#Read in the phytoplankton data
path2phyto <- "D:/"
phyto <- read.table(paste(path2phyto,"a34phy040616.txt",sep=""))


#Convert the class to matrix
phyto <- as.matrix(phyto)


#The data is organised in such a way that 'on land' cells get value -99, R prefers NAs
phyto[which(phyto == -99)][] <- NA


#The orientation is a bit off, so transpose the data and change mirror over latitude (normally you have to look at your data to see how its organised)
phyto <- t(phyto)[,dim(phyto)[1]:1]


#Make a little picture
image(phyto)
#You can see the contours of the UK and the Netherlands & Denmark and even the southern part of Norway in there. Right the area most of our VMS pings are situated as well.


#Although missing from the phyto file, we know what grid cell resolution is associated with this data.
dimnames(phyto) <- rev(list(lat=seq(48.3833,65.9833,1/10)+1/20,lon=seq(-4.1667,13.1667,1/6)[1:83]+1/12))


#The phyto file is currently a matrix, but it would be easier to handle if it would be a data.frame with longitude and latitude columns.
df <- expand.grid(SI_LONG=dimnames(phyto)[[1]],SI_LATI=dimnames(phyto)[[2]])
df$SI_LONG <- anf(df$SI_LONG); df$SI_LATI <- anf(df$SI_LATI)
df$PHYTO <- c(phyto)


#Let's only save those records with data
df <- subset(df,is.na(df$PHYTO)==F)

#Adding at date-time stamp to the phytoplankton data
df$SI_DATIM <- as.POSIXct("25-06-1800 14:00",format="%d-%m-%Y %H:%M")

#Now link the phytoplankton data to the VMS pings
original <- clipObs2Tacsat(tacsat,df,method="grid",control.grid=list(resx=0.1,resy=0.05,gridBbox="obs"),temporalRange=c(-(7*24*60),(7*24*60)),all.t=F)
new      <- clipObs2Tacsat(tacsat,df,method="euclidean",control.euclidean=list(threshold=6),temporalRange=c(-(7*24*60),(7*24*60)),all.t=F)
```

2) In the original code, 21 records were returned, in the new selection, 79 records were returned.
The difference relates to the size of the gridcell specified in the original case. Under those conditions
each gridcell measured approx 7 by 5.5 km. Which meant that VMS pings within that 7 by 5.5km cell were 'linked' to the phytoplankton data.
In the new situation you have selected a much larger area, which is 20km in diameter. And under those conditions, many more VMS
pings could be linked to the phytoplankton data.