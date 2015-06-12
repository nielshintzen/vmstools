# How to link VMS data with other spatial datasets #

## Introduction ##
The number of VMS data requests and analyses have gone up tremendously over the past years. The number of papers written on VMS analyses have followed this increase too. Many of them deal with issues like: how should you clean your data, how should you link VMS and logbook data, what grid size is best to use and also how to calculate your indicators. This already shows that although most of the tools are available, at the advisory level, no consensus has been reached as yet on how exactly to tackle these issues.

At the same time, VMS and logbook data are associated with confidentiality issues, which makes it much more difficult to combine datasets of different nations. Although VMStools presents a framework that can be applied by many different nations at exactly the same manner, also combined with the FishFrame database that might hold VMS data at some point, no such results have been reported yet.

Most of you are here however to do science, and contributing to the discussions outlined above might not be as satisfactory for you. Well, here is a challenge! So far, hardly any papers have discussed the relationship between VMS and covariates as temperature, seabed or wave hight. Why not? Potentially because it is difficult to line all these different datasets up and do the statistics.

Within this practical I'll show you how you can start exploring these interesting aspects, which was one of the core goals to develop VMStools: to quickly deal with the 'un-interesting stuff' of cleaning and linking and get to the explanatory relationships as soon as possible!

## Getting and preparing spatial data ##
One requirement of spatial data is that it contains some sort of position identifier. This could be a shapefile / polygon or a longitude-latitude position(like the start or midpoint of a survey trawl).

In this practical I'll show how to link VMS data with such a shapefile (habitat) and with some phyto plankton data.

### Link VMS with habitat map ###
The good thing about habitat (as in substrate) is that is doesn't change much over time. Once somebody has created a map, you can use it for quite some time and there is no need to update your code. Please now download the [HabitatMaps.zip](http://code.google.com/p/vmstools/downloads/list), store it somewhere on your harddrive (somewhere easy, because you need to link to at later on) and unzip it!

The habitat map is a shapefile and we'll make some plots first to see what it is all about. Thereafter we'll link VMS points up with the habitats and show where most of e.g. the beam trawling has taken place.

```

#Find the path to your habitat map
path <- "D:/HabitatMaps/"
habitat <- readShapePoly(paste(path,"ModelledSeabedHabitats.shp",sep=""))

#We can simply plot it
plot(habitat) #But it doesn't make much sense...
names(habitat)
slotNames(habitat)

#The way the file is constructed is that it has a large number of polygons (167).
#Each of these polygons can be associated with e.g. substrate or salinity.
#The same polygons are used for substrate and salinity but the polygons belonging to the same substrate is not necessarily the same polygons that belong to the same salinity value.

#Let's plot the substrate types around the North Sea
require(maps);require(mapdata)
map("worldHires",xlim=c(-5,10),ylim=c(48,62),fill=T,col="darkgreen");

#What substrate levels do we have?
levels(habitat$substrate) #apparently 7 different substrates

#Define 7 different colors
colrs <- c(
"lightblue","brown","yellow","slategray","green","red",
"blue")

#Turn the polygons into spatial polygons
sp <- SpatialPolygons(habitat@polygons)

#Look at the order of plotting (do I need it: ... not now, but later on yes!)
sp@plotOrder
plot(sp,col=colrs[an(habitat$substrate)],border=colrs[an(habitat$substrate)],add=T)
lcolrs <- length(colrs)+1
legend("bottomleft",legend=c("Land",levels(habitat$substrate)),
lwd=c(1,1,1,1),lty=rep(NA,10),pch=rep(22,lcolrs),ncol=2,
col=rep("black",lcolrs), pt.bg=c("darkgreen",colrs),pt.cex=rep(3,lcolrs),box.lty=1,box.lwd=2,
y.intersp=1.5,title="Substrate type",bg="white",cex=0.5)
map.axes()
```

You might get puzzled from such shapefiles but after a while you'll get a grip on how to handle those! Let's see if we can link VMS pings with substrate type. We must be careful now however. This particular shapefile contains polygons inside polygons. So, if we test if one VMS ping is inside any of the polygons, it might hit 'yes' more than once. So, which one is the right one? Luckily, the file also comes with a plotting order, we need to watch that!

```

#Let's take an example point, e.g. lon,lat = 4.8,54.2. This must be located inside 'Mud to sandy mud'. However, we must be careful here, as the point will be inside more than 1 polygon of the dataset!
#First the easy and automated way
pnt <- c(4.8,54.2)
#Turn the point into a spatial point
spo <- SpatialPoints(coordinates(data.frame(SI_LONG=pnt[1],
SI_LATI=pnt[2])))
#Use the magic 'overlay' function to see in which polygon it is located
idx <- overlay(spo,sp); print(idx)
habitat$substrate[idx]

#Have a look at the map and plot the point, is it really in Mud to sandy mud?
points(pnt[1],pnt[2],pch=19)
```

That wasn't too difficult, and you didn't notice the polygon inside polygon issue at all. That is how nice these things could work. Note that you can have multiple points inside the SpatialPoints. However, if you want / need to do these calculations yourself (because you not only want to perform these analyses, but also want to fully understand what is going on), here is some Do It Yourself code!

```

#Create a vector with zeros representing each of the polygons
# Note that our point can be inside more than one polygon, therefore we have one value (zero or one) for each polygon.
insidePols <- rep(0,length(habitat@polygons))
for(i in 1:length(habitat@polygons)){
for(j in 1:length(habitat@polygons[[i]]@Polygons)){
pol <- coordinates(habitat@polygons[[i]]@Polygons[[j]])
insidePols[i] <- max(insidePols[i],point.in.polygon(pnt[1],pnt[2],pol[,1],pol[,2]))
}
}

#Check inside which polygon our point is situated
idx <- which(insidePols == 1)
print(habitat$substrate[idx])

#This is tricky, we need to define in which substrate the point lies but keep track of the plot order!
habitat$substrate[idx][which.max(order(habitat@plotOrder)[idx])]
```

The main 'engine' of most of these calculations is the point.in.polygon function which is very powerful, but the way it iss coded above is very inefficient for just one point!

Time to find some relationships. Let's take some beamtrawl fishers only and see what type of substrate they normally trawl. Do you have a gut-feeling already?

```

#Load the data, merge it, sort and filter and make a subset
data(tacsat);data(eflalo)
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]
tacsat  <- subset(tacsatp,LE_GEAR == "TBB")
tacsat  <- sortTacsat(tacsat)
tacsat  <- filterTacsat(tacsat,c(2,6),hd=NULL,remDup=T)

#All the preparation is done now, let's overlay these VMS points with the habitat map (I am interested in substrate)
spo     <- SpatialPoints(coordinates(tacsat[,c("SI_LONG","SI_LATI")]))
sp      <- SpatialPolygons(habitat@polygons)
idx     <- overlay(spo,sp) #Get a cup of coffee, this might take a while...!
tacsat$SUBSTRATE <- habitat$substrate[idx]

table(tacsat$SUBSTRATE)
```

### Link VMS with EEZ map ###
Because the habitat map is so big, and might cause memory problems, we also do very similar things with an EEZ shapefile which is much smaller in size.

Download the EEZ shapefile [EEZ.shp](http://code.google.com/p/vmstools/downloads) and place it in a folder you can find easily and unzip the file. We will read the file in, make some simple plots of the total shapefile and extract the Dutch EEZ out of it. Thereafter we will see and plot which tacsat points are located inside the EEZ and which ones outside. A practical example could be that I am interested to study which fraction of the total Dutch fleet effort is allocated inside and outside the EEZ.

```

#Set path to EEZ file
path <- "C:/tmpTacsat/"

#Read all EEZ of the world
EEZ         <- readShapePoly(paste(path,"World_EEZ_v2_20081031.shp",sep=""));

#Plot the EEZ's
plot(EEZ)

#Take the EEZ of the Netherlands
which(EEZ$Country == "Netherlands") #165

#To substract the Dutch EEZ out of the full file, I first need to turn the EEZ into a SpatialPolygons
EEZ <- as(EEZ,"SpatialPolygons")

#Now select EEZ polygon 165
EEZNLD <- EEZ[165]

plot(EEZNLD)

#Overlay the EEZ with the tacsat points, I first need to convert my tacsat points to SpatialPoints however.
data(tacsat)

sptacsat <- SpatialPoints(coordinates(tacsat[,c("SI_LONG","SI_LATI")]))

idx      <- overlay(sptacsat,EEZNLD)

#Plot the Dutch EEZ once more, and add the points inside the EEZ in blue and outside in red.
plot(EEZNLD)
points(tacsat$SI_LONG[which(idx==1)],tacsat$SI_LATI[which(idx==1)],col="blue")
points(tacsat$SI_LONG[which(is.na(idx)==T)],tacsat$SI_LATI[which(is.na(idx)==T)],col="red")

map.axes()
```

### Link VMS with phyto plankton map ###
Some datasets however, we'd like to link to VMS, change considerably over time. Think for example about seabirds following fishing vessels, the benthic community composition during a year or phytoplankton abundance in the North Sea. This requires us to 'overlay' the VMS and other datasets in both time and space. This is obviously different from the exercise above where we were only overlaying space.

As an example, we will use the phytoplankton distribution data that can be downloaded as [a34phy040616.txt](http://code.google.com/p/vmstools/downloads/list). Again, download it to a place from where you can easily read it in. This dataset only represents 2 weeks of data, and can therefore only be associated with 2 weeks of VMS pings too! Let's get started!

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
```

The data needs to get a time stamp too, because it only represents 2 weeks (and allows you to add more weeks of data below the current data in the same dataframe). In that way, you can link your VMS with e.g. a whole year of phytoplankton data at once.

```

#Adding at date-time stamp to the phytoplankton data
df$SI_DATIM <- as.POSIXct("25-06-1800 14:00",format="%d-%m-%Y %H:%M")

#Now link the phytoplankton data to the VMS pings
res <- clipObs2Tacsat(tacsat,df,method="grid",control.grid=list(resx=0.1,resy=0.05,gridBbox="obs"),temporalRange=c(-(7*24*60),(7*24*60)),all.t=F)
```

That was just linking the two datasets, however, we would like do some sort of analyses with it. Perhaps associate each VMS ping with a phytoplankton value! Keep well in mind however that we do not link phytoplankton to VMS but VMS to phytoplankton! So only those VMS records that are within the 2 weeks of our phytoplankton data, and appear at the right spot will be linked to our phytoplankton data. In statistics, a model will look like this: phytoplankton ~ VMS + .... and not VMS ~ phytoplankton!

```

head(res[[1]]) #The total observation dataset, now gridded with grid index value
head(res[[2]]) #All the linked VMS points, also with grid index value
tacsatLinked <- res[[2]]
tacsatLinked$PHYTO <- res[[1]]$PHYTO[match(res[[2]]$GR_ID,tacsatLinked$GR_ID)]

#A completely arbitrary picture of the results
# Do they fish harder when there is more phytoplankton?
plot(tacsatLinked$SI_SP~tacsatLinked$PHYTO,xlab="Phytoplankton abundance",ylab="Speed (knots)",pch=19)

#And which phytoplankton datapoints got linked up in the end:
plotTools(tacsatLinked,level="gridcell",xlim=c(-4,9),ylim=c(51,60),gridcell=c(0.1,0.05),control.tacsat=list(clm="PHYTO"),col=heat.colors(9))
```

### Exercise 1 ###
  1. There is an option to change from grid to a 'distance' based calculation in the 'clipObs2Tacsat' function. Perform the same analyses as above, but now with method="euclidean" and threshold = 10km.
  1. How many records are returned in the first analyses and how many in the second, can you explain the difference