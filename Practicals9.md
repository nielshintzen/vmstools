# Pictures #

## Introduction ##
In many cases, VMS and Logbook data are used to create maps of fishing effort in specific areas or over time. Although R facilitates making plots,
making one that looks good may be hard work. Here we present one potential way of doing it, and you can tailor the specifications to your own
needs. See it as a template which you can use over and over, that is at least how it's currently being used.

## Picture settings ##
Setting your picture settings in advance, at the top of your script, might seem a lot of extra work.
However, once you've spefied all your settings at once, you'll make sure that in each figure you design
afterwards you will use the exact same settings.

```

library(vmstools)
library(maps)
library(mapdata)
library(RColorBrewer)

cl          <- 1.1  #cex.lab
ca          <- 1    #cex.axis
fonts       <- 2    #font
xl          <- list(label="Longitude",font=fonts,cex=cl) #x-label
yl          <- list(label="Latitude",font=fonts,cex=cl)  #y-label
zl          <- list(font=fonts,cex=cl) #z-label (if any)
colintens   <- brewer.pal(6,"YlOrRd")  #colour for fishing intensity
colland     <- brewer.pal(9,"PiYG")[8] #colour of land
colgrey     <- brewer.pal(9,"Greys")   #colour of grey shades
figtype     <- "tiff"                  #figure extension type
parmar      <- rep(2,4)                #panel settings
paroma      <- (c(6,6,2,2)+0.1)        #panel settings
reso        <- 1                       #dpi setting (1=100dpi)
```

### Specifying the areas to plot ###
The next step is setting the spatial extend of the figure, creating a grid to project the
VMS data to and potentially loading a few additional area sets (such as N2000 areas)
The code below shows you how to create such a grid, how to convert your VMS data in a way that it
can easily be projected on a grid and it defines some other areas of interest.

```

#- Define grid cell area
resx        <- 0.25
resy        <- 0.125

#-Load the data and convert your data into a 'SpatialPointsDataFrame'
data(tacsat); data(eflalo)
coords      <- SpatialPointsDataFrame(cbind(x=an(ac(tacsat$SI_LONG)),y=an(ac(tacsat$SI_LATI))),data=tacsat)
```
A 'SpatialPointsDataFrame' is the same as a 'SpatialPoints' object, only this one allows you to store per coordinate a lot of data inside a data.frame.

```

data(ICESareas)
data(europa)

#-Define the area of interest and some reference area
areaInt    <- ICESareas[grep("IV",names(ICESareas))[1:3]]
areaRef    <- ICESareas[grep("IVb",names(ICESareas))[1]]

#-Obtain outer region of my areas and VMS positions. This helps to create maps later on with the same dimensions.
bbox        <- cbind(bbox(areaInt),bbox(areaRef))
spatBound   <- list(xrange = c(floor(range(bbox["x",])[1]),ceiling(range(bbox["x",])[2])),
yrange = c(floor(range(bbox["y",])[1]),ceiling(range(bbox["y",])[2])))
grd         <- createGrid(spatBound$x,spatBound$y,resx,resy,type="SpatialGridDataFrame",exactBorder=T)
```
We have loaded a file called 'europa' before, below we convert it into a SpatialPolygons object which allows easy plotting
```

#This is a bit slow though
eurPols     <- lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),
function(x){data.frame(SI_LONG=subset(europa,SID==x)$X,SI_LATI=subset(europa,SID==x)$Y)}))
```

All areas are defined now. However, once we start plotting the results, quite often
we find out that we need some colour scale to present the results at. If for example
the data spans from 0 to 1000, some steps like 0-10,10-25,25-50,50-100,100-200,200-400 might
look good. Below we define these ranges first, which makes our lives easy when plotting and
adding legends to these plots.

```

cutbreaksval  <- list(ALL = c(-1,0,10,25,50,100,200,400))
legval        <- list(ALL = c("0","0 <= 10","10 <= 25", "25 <= 50","50 <= 100","100 <= 200","200 <= 400"))
#- Potentially, divide your data by a certain constant and add this constant to the legend title
valdiv        <- 1000 * 60 #combination of converting from minutes to hours and getting a legend value per 1000 hours
unitval       <- c('x 1000 hours per year')
```

## Creating the first figure ##
Since we've defined all figure settings and areas, we can start to make some simple plots.
Let's start with a figure of the study areas

```

tiff(filename=paste("D:/output/","Map_interestAreas.tiff",sep=""),compression = c("none"),bg="white",width = (17/2.54), height = (25/2.54), units = "in", pointsize = 12,res=reso*100)

plot(1,1,col="white",xlim=spatBound$xrange,ylim=spatBound$yrange,xlab="",ylab="",las=1,cex.lab=xl$cex,font=xl$font,
asp=1/lonLatRatio(mean(spatBound$xrange),mean(spatBound$yrange)))
plot(areaInt,add=T,col="yellow")
plot(areaRef,add=T,col="lightblue")
plot(eurPols,add=T,col=colland)
axis(1);axis(2,las=1); box()

mtext(xl$label,side=1,outer=T,line=-3,at=0.5,font=xl$font,cex=xl$cex)
mtext(yl$label,side=2,outer=T,line=-1.5,at=0.5,font=yl$font,cex=yl$cex)
dev.off()
```

Note in the picture above that ICES area IVb was first plotted in yellow and thereafter, on top of that,
in blue.

Since we've now got the study area available, we can add some simple figures on fishing intensity.
In principle, the code to plot swept area, value of an area, trawling frequency or intensity is all the same.
The only difference will pop up in the column you'll try to aggregate / summarize over.

```

tiff(filename=paste("D:/output/","valueMap_allGearsLarge.tiff",sep=""),compression = c("none"),bg="white",width = (17/2.54), height = (25/2.54), units = "in", pointsize = 12,res=reso*100)

#- Reset values
grd@data[] <- 0

#-Create column to aggregate over (most often, this column already exists and is output from your previous analyses)
tacsat                        <- intervalTacsat(tacsat,level="vessel",fill.na=T)

idx                           <- over(as(coords,"SpatialPoints"),as(grd,"SpatialGrid"))
tacsat$gridID                 <- idx

#- Here we aggregate data to the grid cell. You can aggregate any column you like, we use INTV as an example here.
grd@data[names(table(idx)),1] <- aggregate(tacsat$INTV,by=list(tacsat$gridID),FUN=sum,na.rm=T)$x


plot(1,1,col="white",xlim=spatBound$xrange,ylim=spatBound$yrange,xlab="",ylab="",las=1,cex.lab=xl$cex,font=xl$font,
asp=1/lonLatRatio(mean(spatBound$xrange),mean(spatBound$yrange)))
coordGrd                      <- coordinates(grd)[an(names(table(idx))),]
#-Here we turn each of the grid cells into a polygon. All polygons together make the picture
grdPols                       <- lonLat2SpatialPolygons(lst=lapply(as.list(1:nrow(coordGrd)),
function(x){data.frame(SI_LONG=c(coordGrd[x,"s1"]-resx/2,rep(coordGrd[x,"s1"]+resx/2,2),coordGrd[x,"s1"]-resx/2),
SI_LATI=c(rep(coordGrd[x,"s2"]-resy/2,2),rep(coordGrd[x,"s2"]+resy/2,2)))}))
```
We have specified the binns we use for the legend and colours above, but it might be handy if you see where I got the range in values from.
Below, you can see the 'grd@data[an(names(table(idx))),1]' piece of code. This selects only those coordinates from the data.frame that are
associated with a value > 0. All in all, this is a simple vector and allows me to see what values it spans.
```

#Look at value ranges:
range(grd@data[an(names(table(idx))),1])

#-Here we cut the range of values up into colour binns
cols                          <- c("white",colintens)[cut(grd@data[an(names(table(idx))),1]/valdiv,breaks=cutbreaksval$ALL)]
plot(grdPols,col=cols,add=T,border='transparent')
plot(eurPols,add=T,col=colland)
axis(1);axis(2,las=1); box()

#-Add a legend
legend(x='topright',fill=c('white',colintens),legend=legval$ALL,bg='white',title=unitval,box.lty=1)

#- Add axis and title
title(main="Fishing itensity example per 1/16th ICES square",outer=F,cex=ct)
mtext(xl$label,side=1,outer=T,line=-3,at=0.5,font=xl$font,cex=xl$cex)
mtext(yl$label,side=2,outer=T,line=-1.5,at=0.5,font=yl$font,cex=yl$cex)
dev.off()
```

## Shape files ##
If you prefer not to plot it straight into R, but rather export it to e.g. ArcGIS, then below
you'll find a routine to export the 'grdPols' we created above to a shapefile.

```

library(rgdal)
grdPolsDF              <- as(grdPols,"SpatialPolygonsDataFrame")
grdPolsDF@data         <- data.frame(value=grd@data[an(names(table(idx))),1],color=cols)
proj4string(grdPolsDF) <- CRS("+proj=longlat +ellps=WGS84")
dir.create("D:/output/shp_files/")
setwd(file.path("D:/output/","shp_files"))
writeOGR(grdPolsDF, dsn = '.', layer = 'Fishing_intensity_tacsat', driver = "ESRI Shapefile")
```

### Exercise 1 ###
  * Try to open the shapefile in your GIS package and investigate if it looks the same.

# Tables #

## Introduction ##
Creating tables out of the VMS and Logbook data is often based on 'point in polygon' types of
calculations where GPS positions inside a refenence area or area of interest are aggregated.
Here we show how you could also merge logbook records with VMS records in the same type of
analyses.

Let's start with calculating the time spend in areas IV and IVb in specific by VMS and Logbook data.
```

data(tacsat); data(eflalo)

#Calculate effort associated with each VMS ping
tacsat          <- intervalTacsat(tacsat,level="vessel",fill.na=T)
tacsat          <- tacsat[which(is.na(tacsat$SI_LONG)==F & is.na(tacsat$SI_LATI)==F),]
#Calculate the effort in an eflalo dataset
eflalo          <- effort(eflalo,unit="mins")
#Because some trip numbers are replicated, we take only unique trip numbers
# (a single trip can be associated with multiple eflalo records, that means that also time out of harbour and in harbour are replicated)
eflalo          <- eflalo[!duplicated(eflalo$FT_REF),]
#Assign GPS positions to eflalo
eflalo$SI_LONG  <- ICESrectangle2LonLat(eflalo$LE_RECT)$SI_LONG
eflalo$SI_LATI  <- ICESrectangle2LonLat(eflalo$LE_RECT)$SI_LATI
eflalo          <- eflalo[which(is.na(eflalo$SI_LONG)==F & is.na(eflalo$SI_LATI) == F),]

coordsTacsat    <- SpatialPoints(tacsat[,c("SI_LONG","SI_LATI")])
coordsEflalo    <- SpatialPoints(eflalo[,c("SI_LONG","SI_LATI")])

#Determine tacsat GPS position inside area of interest (IVa, IVb and IVc)
# and determine which midpoint of eflalo rectangles are within the area of interest
inTacsat        <- over(coordsTacsat,areaInt)
inEflalo        <- over(coordsEflalo,areaInt)

#Calculate here the total minutes spend by vessels in these areas
sum(tacsat$INTV[      which(is.na(inTacsat)==F)],na.rm=T) +
sum(eflalo$LE_EFF_LOG[which(is.na(inEflalo)==F)],na.rm=T)
```

Easier would it be if we could combine tacsat and eflalo straight away. This especially
pays of if you start to divide ICES rectangles up into e.g. 1/4th 1/16th or even small
areas. The midpoints of each of these 1/4th or 1/16th ICES rectangles can be treated
similarly.

```

eflaloTrim      <- eflalo[,c("SI_LONG","SI_LATI","LE_EFF_LOG")]
colnames(eflaloTrim)[3] <- "INTV"
comb            <- rbind(tacsat[,c("SI_LONG","SI_LATI","INTV")],
eflaloTrim)
coordsComb      <- SpatialPoints(comb[,c("SI_LONG","SI_LATI")])

inComb          <- over(coordsComb,areaInt)

#Calculate the total minutes spend by vessels from tacsat and eflalo combined
sum(comb$INTV[which(is.na(inComb)==F)],na.rm=T)

#Now do it for IVb only
inComb          <- over(coordsComb,areaRef)
sum(comb$INTV[which(is.na(inComb)==F)],na.rm=T)

#Or we define an arbitrary polygon ourselves
arbPol          <- lonLat2SpatialPolygons(SI_LONG=c(2,3,3.5,2),SI_LATI=c(53,54,56.4,53))
plot(1,1,col="white",xlim=spatBound$xrange,ylim=spatBound$yrange,xlab="",ylab="",las=1,cex.lab=xl$cex,font=xl$font,
asp=1/lonLatRatio(mean(spatBound$xrange),mean(spatBound$yrange)))
plot(eurPols,add=T,col=colland)
plot(arbPol,col=2,add=T)
axis(1);axis(2,las=1); box()

inArbPol        <- over(coordsComb,arbPol)
sum(comb$INTV[which(is.na(inArbPol)==F)],na.rm=T)