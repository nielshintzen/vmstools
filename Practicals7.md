# Practical 7: Output to Fishframe and GIS #

## Output to FishFrame ##

  * Here we describe a method due to Bastardie et at. (2010) for simultaneously coupling VMS and logbook data and then distributing the landings or their cash values among VMS positions estimated to be actively fishing.

  * The first step is to work out the proximity of each VMS data point to each fishing port.  Linking each VMS trip to each logbook trip is done by searching for the nearest trip among all logbook trips based on the temporal mid-point of trips.  Note: that the method described here differs from that due to Hintzen and described in Practical 3.

  * After linking VMS and logbook data, gear types can be assigned to each VMS position enabling fishing activity, or otherwise, to be identified using the segmented regression approach described earlier.

  * The next step is to evaluate the consistency between each of the two separate data sources (VMS vs logbook). This process is described in detail by Bastardie et al. (2010).

  * Finally the landings (in both weight and cash value) declared in the logbooks are allocated to each VMS position.

  * The data can then be reshaped and aggregated into any other format and the example used here is FishFrame. It is worth noting that R-scripts like the one described below for our 'Atlantis' data can of course be run on VMS and logbook data for any different country to produce aggregated datasets in common formats that can be combined and exchanged without violating data protection laws.

  * The scripting method ensures that each dataset is processed in an identical way and is repeatable.


## Getting started ##
As usual we have to add the inbuilt data sets to R's search path, install the vmstools library, and 'clean' the data prior to any analysis.

Loading the vmstools library:

```

library(vmstools)
```

Hereafter, as we've seen before, we need to load the logbook, VMS and harbour data to have access to these data sources.

```

#Load the data
data(eflalo)
data(tacsat)
data(euharbours)
```

As usual, and especially necessary for your own data, format the data:

```

#Format and sort the data
eflalo <- formatEflalo(eflalo)
tacsat <- formatTacsat(tacsat)
tacsat <- sortTacsat(tacsat)
```

The code that turns the data into FishFrame data has a slightly different setup than the other functions you've seen so far.  The manner via which logbook and VMS data are linked is different and it requires a pre-definition of points in harbour. Here we test whether VMS positions are in harbour or not:

```

#Find out in which harbour they are located
tacsat$SI_HARB <- NA
harbours$Description <- harbours$harbour
tacsat$SI_HARB <- pointInHarbour(lon=anf(tacsat$SI_LONG),
lat=anf(tacsat$SI_LATI),
harbours=harbours,
rowSize=30, returnNames=TRUE,saveHarbourList=F)

#Create a separate object with just the harbour information:
inHarb <- tacsat$SI_HARB
inHarb <- replace(inHarb, !is.na(inHarb), 1)
inHarb <- replace(inHarb, is.na(inHarb), 0)
inHarb <- as.numeric(inHarb)
head(inHarb) #What did just happen?
```

Once we've defined which points are in harbour, the functionality requires some sort of trip identifier:

```

tacsat$SI_FT <- 1 # init
idx <- which(inHarb==0)
tacsat[idx,"SI_FT"] <- cumsum(inHarb) [idx] # add a SI_FT index
```

Linking logbooks and VMS is defined here based on mid-points of trips and therefore requires points outside of a harbour but also the departure and arrival dates.

```

#A way to get a trip identifier to the tacsat dataset, some clever coding, but perhaps a bit more difficult to understand. Have a look at what the 'diff' function really does
startTrip <- c(diff(tacsat[,"SI_FT"]), 0)
endTrip   <- c(0, diff(tacsat[,"SI_FT"]))
tacsat[which(startTrip>0),"SI_FT"]  <-  tacsat[which(startTrip>0)+1,"SI_FT"]
tacsat[which(endTrip<0),"SI_FT"]    <-  tacsat[which(endTrip<0)-1,"SI_FT"]
tacsat <- tacsat[which(inHarb==0 |  startTrip>0 |  endTrip<0),]
```

Now assign activity state (e.g. fishing v steaming) to each ping. Note that this is optional. If detectFishing = TRUE in mergeEflalo2Pings then the function/code will do this part automatically. You might want, for example, to apply your own rules to identify fishing. You can do it here by setting detectFishing = FALSE.

```

#Add activity to the dataset
tacsat$SI_STATE <- 2 # init (1 = fishing; 2 = steaming)
tacsat$SI_STATE [(tacsat$SI_SP>4 & tacsat$SI_SP<8)] <-1
```

Sometimes you get naming mistakes in eflalo and below is some example code that might be useful to deal with these:

```

#    if(!match('LE_MET_level6',colnames(eflalo))>0){
#      eflalo$LE_MET_level6 <- eflalo$LE_MET
#    }
```

In this example we run the code for a given set (50 random ones) of vessels in the year 1800. Since detectFishing = TRUE (see below) then the code is also detecting fishing activity using the built-in function segmentTacsatSpeed(). [that this will overwrite the existing SI\_STATE](Note.md). First we set up some variables we need,ie. a year variable, and a vector for the unique boat names and also the total number of boats that we intend to output data for.

```

#get a year column in the data
tacsat$year <- as.numeric(format(tacsat$SI_DATIM, format = "%Y"))

#Making a selection of boats
how.many.boats.in.tacsat <- length(unique(tacsat$VE_REF     [tacsat$year==1800]))  #Take data for 1800
no.of.boats              <- 50 # lets's do 50 boats
which.boats <- unique(tacsat$VE_REF)[seq(1,how.many.boats.in.tacsat,length=no.of.boats)]
```

This code writes data to the c:\\ drive so it's a good idea to clean it out first to avoid confusion:

```

unlink("c:\\output\\*") #delete everything in c:\\output to avoid confusion.
```

Now run the main program that does all the work. You will note pictures popping up etc. These can all be found in c:\\output and describe the matching of the VMS to logbook data and how fishing was identified for each vessel.

```

mergeEflalo2Pings (eflalo=eflalo, tacsat=tacsat, vessels=which.boats,
general=list(output.path=file.path("C:","output"),
visual.check=TRUE, detectFishing=FALSE, speed="segment",
what.speed="calculated"))

```

Another option might be to run the code for ALL the vessels you have available. The code to do this is described below should you want to try...

```

#>mergeEflalo2Pings (eflalo=eflalo, tacsat=tacsat,
#                     general=list(output.path=file.path("C:","output"),
#                     visual.check=TRUE, detectFishing=TRUE, speed="segment",
#                     what.speed="calculated"))
```

Tidy up the results in your R directory:

```

gc(reset=TRUE)
```


Similarly you might want to load the merged output table for one vessel, say vessel "10" for the year 1800 and check whether all the landings in eflalo (logbook) have actually been allocated to VMS positions, which we do below:

```

#Load some of the results we got above
load(file.path("C:","output","merged_10_1800.RData"))
merged.10 <- merged

#Sum over the plaice landings
sum(tapply(anf(merged$LE_KG_PLE), merged$flag, sum, na.rm=TRUE))
sum(eflalo[eflalo$VE_REF=="10","LE_KG_PLE"], na.rm=TRUE) #These numbers should be the same.
```

Alternatively we can bind up all 50 vessels keeping the file sizes manageable by selecting just a few species:

```

bindAllMergedTables (vessels=which.boats, a.year = "1800",
species.to.keep=c("PLE","COD","SOL","DAB","TUR"),
folder = file.path("C:","output"),
all.in.one.table=TRUE)
```

We then load up the merged output table for all 50 vessels:

```

load(file.path("C:","output","all_merged__1800.RData"))
```

Switch graphics off

```

graphics.off()
```

Time to map the data!

```

# Select columns of interest for convenience.

df1<- all.merged[, c("SI_LATI","SI_LONG","LE_KG_COD","LE_KG_PLE","LE_KG_SOL","LE_KG_DAB","LE_KG_TUR")]

#Need to turn my longitude and latitude into factors here, this is tricky, and we do it for a very specific reason (another function that needs it, so don't try this at home (unless you want to copy what we do here :-))

df1$SI_LONG <- anf(df1$SI_LONG)
df1$SI_LATI <- anf(df1$SI_LATI)
df1 <-   df1[ !is.na(df1$SI_LATI),];
df1 <-   df1[ !is.na(df1$SI_LONG),]

#Chose species to look at
what.species.to.map <- "LE_KG_PLE"

df2<-df1[-which(is.na(df1[,what.species.to.map])),] # df2 are the data for plaice only

# Note: we use vmsGridCreate here which is being phased out and replaced by createGrid. plotMap = TRUE
# will send the data to the R-graphics windows.

vmsGridCreate(df2,nameLon="SI_LONG",nameLat="SI_LATI", we = 3, ea = 6, so = 50, no = 54,
nameVarToSum = what.species.to.map,cellsizeX =0.05,plotMap=TRUE,
cellsizeY =0.05,  legendtitle = paste(substr(what.species.to.map,7,9)," landings (kg)")
, plotPoints =TRUE, breaks0=c(1,4,8,16,32,64,100000))

```

Now convert to FishFrame format which might take some time running. By default this function will keep all the species in the output table

```

# First bind up all the tables in c:\\output
tmp <- bindAllMergedTables (vessels=which.boats,
a.year=1800,
species.to.keep=character(),
folder = file.path("C:","output"),
all.in.one.table=FALSE)
#Now create FishFrame data. Have a look at the Downloads section to get an explanation of fishframe format and have a look at the file 'ff'.

ff  <- pings2Fishframe (general=list(output.path=file.path("C:","output"),
a.year=1800, a.country="NLD",degree=0.05) )
```

Note that the file of interest (FishFrame format) will be in c:\\output and is called for this run of the script, NLD1800\_VD.csv. Open it in Wordpad or similar and have a look at it

## Output to more formal GIS: ##

Add in the arguments for a grid below:

```

# Set up desired grid size and range

grid.size <- 0.06
we <- -3
ea <- 10
so <- 49
no <- 59

#Directory to put output file in and the name of the file
setwd('c://output')
outfile <- 'fishlandings.asc'

#Run the function vmsGridGreate, again, which has an option to output an ASCII grid that is in ArcView format but can also be read by other GIS programs. Set plotMap = F now.

output.mat<-vmsGridCreate(df2, we=we,ea=ea,so=so,no=no,nameLon = "SI_LONG", nameLat = "SI_LATI",plotMap=F,
nameVarToSum="LE_KG_PLE",cellsizeX = grid.size, cellsizeY = grid.size,
legendx="bottomright",numCats=9,plotPoints=F,legendncol=2,legendtitle="kgs",outGridFile=outfile)
```

A free and actually very good GIS system is Google Earth which you can also use easily with our vmstools library (note that GE avoids the problems of 'projections' because it is all actually done on the surface of a sphere). If Google Earth is not already installed on your machine, do so now!

```

library(nlme);library(spatstat) # Install them if you have not already.
# Create an output filename
fn<-'landings'
what.quantity = 'kgs' # What quantity are we plotting ?

# Re-organise the output from vmsGridCreate:
dd <- output.mat@grid;
d1 <- dd@cells.dim[1];
d2 <- dd@cells.dim[2];
fishing <- output.mat@data$fishing;
mat <- matrix((fishing),byrow=TRUE,ncol=d1,nrow=d2);
mat <- t(mat[d2:1,]);
log.it <- T; # do you want to log transform your data before plotting them ?
bbox <- output.mat@bbox; # Extract the bounding box
xxx <- seq(bbox[1,1],bbox[1,2],length=d1); # Create a vector of longitudes for the 'image' function.
yyy <- seq(bbox[2,1],bbox[2,2],length=d2); # Create a vector of latitudes for the 'image' function.
rr <- range(mat[mat!='-Inf'],na.rm=T);     # Get the range of the data.
labs<-seq(rr[1],rr[2],length=9);
image(xxx,yyy,mat,zlim=c(rr[1],rr[2]),xlab="",ylab="",col=rainbow(9)); # Plot to R to check things are ok.

# Format the data for plotting in Google Earth.
gd <- list(x=xxx,y=yyy,z=mat);
gd.1 <- as.SpatialGridDataFrame.im(as.im(gd));
# create right format for GE_SpatialGrid
proj4string(gd.1) <- CRS("+proj=longlat +datum=WGS84");
vms.kml <- GE_SpatialGrid(gd.1);
#tf <- tempfile(tmpdir=getwd())

imagefile <- 'kml.png'
# create name for png file.
kmlfile   <-  paste(fn,'.kml',sep="")

# Create a png file to load into Google Earth.
png(file=imagefile, width=vms.kml$width, height=vms.kml$height, bg="transparent",res=576);
par(mar=c(0,0,0,0), xaxs="i", yaxs="i",cex=.25);
image(as.image.SpatialGridDataFrame(gd.1[1]), col=heat.colors(9),xlim=vms.kml$xlim, ylim=vms.kml$ylim);
kmlOverlay(vms.kml, kmlfile=kmlfile, imagefile=imagefile, name=what.quantity);
# legend(x='bottomright', legend=as.character(round(labs)), pch = 22,
# border =heat.colors(length(labs)), pt.bg=heat.colors(length(labs)),
# title=what.quantity, ncol=2,bg="white",pt.cex=1 ) # This is not working properly. Perhaps 'color.legend' might work better?
dev.off()

# There is also a function that does this called 'Grid2KML' which is in the old part of the repository (https://vmstools.googlecode.com/svn/trunk/vmstools/R/old/Grid2KML.r). It needs updating but you can try it if you like.

Grid2KML(output.mat=output.mat,what.quantity = 'kgs',kmlfile=paste(fn,'.kml',sep=''),imagefile=paste(fn,'.png',sep=''))

```

Now if you go to the directory, c:\\output, and doubleclick on the kml file there Google Earth will open automatically upon which your data will be superimposed!


### Exercises ###
  1. Go and install QGIS (http://www.qgis.org/) which is free and useful GIS system.
  1. Find an .asc file in c:\\output output by vmsGridCreate and load into QGIS as a raster layer.
  1. In QGIS apply a color scale to the monochrome grid that you should be able to see.
  1. In QGIS convert this 'raster' file to a 'vector' file format (e.g. shapefile), save it to disk and then import into R and plot.

Note that the idea here is just to show that vmstools need not be used in isolation from more formal GIS, and that many spatial file formats can be handled in R. Furthermore you don't have to pay for good GIS systems either! Udig is another good option by the way http://udig.refractions.net/.















