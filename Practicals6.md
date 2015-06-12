# Practical 6: Calculating indicators #

## Introduction ##
All of the main marine policy frameworks in Europe have adopted an ecosystem approach to  managing their anthropogenic activities and conserving their natural resources. For effective ecosystem-based management (EBM), indicators are required that describe the pressures affecting the ecosystem, the state of the ecosystem, and the response of managers. To implement EBM successfully it is also critical that we understand the links between these indicators, so that the policy objectives on the state of the ecosystem can be achieved by managing human pressures through an appropriate response.

The EU Data Collection Framework (DCF) standardizes three indicators to analyse fishing activity. They are summarised as follow:

  * DCF Indicator 5: Distribution of fishing activities. The spatial extent of fishing activity based on the total area of grids within which VMS records were obtained, each month;

  * DCF Indicator 6: Aggregation of fishing activities. The extent to which fishing activity is aggregated based on the total area of grids within which 90% of VMS records were obtained, each month.

  * DCF Indicator 7: Areas not impacted by mobile bottom gears. The area of seabed that has not been impacted by mobile bottom fishing gears in the last year. Could be reported annually and would state the total proportion of the area by depth strata in each marine region.

Although the definitions of these indicators have been around for some time already, and tools are even developed to calculate them, there is still room for alternative approaches.

Within this practical we will guide you through the simplest ways to calculate the indicators, and a home-made approach which allows for more flexibility.

### Some background on the indicators and their calculation ###
These indicators aggregate the tacsat point data into a gridded data frame using the functions mapGrid and vmsGridCreate, and therefore the resolution of the grid (cell size) must first be defined.

**DCF 5** calculates the total area of a grid of cells with fishing activity which is above a minimum threshold of number of pings or number of fishing hours (if the tacsat data contain a field with time interval between two points called SI\_INTV, then the threshold will be a minimal number of hours, otherwise it will be a minimal number of points). The area of each cell is calculated with the function 'surface', either via a fast and rough method using a trapezoid  approximation (option "Trapezoid"), or via a more precise but slow method using a Universal Transverse Mercator projection (option "UTM"). The first method is fine for big grids of small cell sizes, the second one is better for large cell sizes. This total fishing area is processed by month.

**DCF 6** also calculates the total area of a grid with fishing activity but keeps only 90% of the points by discarding the outer 10% of points (or any other specified percentage). It uses the function 'tacsatMCP' adapted from the 'aspace' library. This function draws a minimum convex polygon around the central points to keep. These points are then gridded and the total area of the cells is calculated with the 'surface' function with the same optional methods as DCF 5. The total fishing area is processed by month.

**DCF 7** calculates the total area of a specified polygon not impacted by mobile bottom gear. The tacsat data must, therefore, have first been merged with the logbooks in order to have a gear code (or others) for each VMS point. The indicator needs a list of gear codes to include as mobile bottom gears (if empty, all the points will be included). The specified area to be processed is a 'polygon shapefile'. This polygon (or group of polygons) is then clipped with the Europe polygon to be sure that the indicator won't include land in its area calculation. If no shapefile is defined, the area of the bounding box containing all the VMS points will be considered. The result is the area of the polygon, less the area of the grid where fishing activity occurs. The VMS pings are gridded with an optional threshold in either minimum of fishing hours or minimum number of points (see DCF 5). The area of each grid cell is calculated with the 'surface' function (see DCF 5 or DCF 6).


## Calculating indicators the easy way ##
Again, we have to start with the tacsat and eflalo data, make sure they are free from the easy-to-spot errors and link them up. As we do this for example purposes only, we skip the cleaning steps and move to linking the tacsat and eflalo datasets immediately.

```

#Load the data, merge the datasets together and assign gear to tacsatp
data(tacsat); data(eflalo)
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]

#In general, most of these indicators only make sense for bottom gears, so let's make a selection
gear <- c("OTB","PTB","TBB","OTM")
tacsatp <- subset(tacsatp,LE_GEAR %in% gear)

#Again, quick and dirty filtering activity
tacsatp <- filterTacsat(tacsatp,st = c(2, 6), hd = NULL, remDup = T)
tacsatp <- sortTacsat(tacsatp)

#Calculate the interval rate between pings
tacsatp <- intervalTacsat(tacsatp,level="trip",fill.na=TRUE)
tacsatp$SI_INTV <- tacsatp$INTV

#Define an area of interest (prevent having an influence of VMS pings in a completely different area)
spatBound <- list(x=c(-5,10),y=c(48,62)) #Spatial boundaries of approx the North Sea

#Just some code to select only those tacsat records which are within these boundaries
idxlon              <- which(tacsatp$SI_LONG >= spatBound$x[1] & tacsatp$SI_LONG <= spatBound$x[2])
idxlat              <- which(tacsatp$SI_LATI >= spatBound$y[1] & tacsatp$SI_LATI <= spatBound$y[2])
tacsatp             <- tacsatp[idxlon[which(idxlon %in% idxlat)],]

#Let's make a quick graph
plotTools(tacsatp,level="ICESrectangle",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,log=FALSE,gridcell=c(0.1,0.05),color=NULL,control.tacsat=list(clm="SI_INTV"))
```

Ok, all the data has been cleaned, the interval rate has been calculated (optional, otherwise VMS pings will be taken as a measure) and roughly the North Sea area has been selected. Now onto calculating the indicators!

```

#Set a working directory where the files will be written to disk
setwd("C:/tmpTacsat/")
ind5                <- indicators( indicatorNum=5,
tacsatp,
minThreshold=0, # percentage of points to be included in the Minimal Convex Polygon
cellresX=0.1,
cellresY=0.05,
calcAreaMethod="Trapezoid",
plotMapTF=T,
exportTableName="indicator5_01005",
exportGridName="indicator5_01005")
print(ind5)

#Let's change the resolution of the calculated grid and see how the indicator values change
ind5_2                <- indicators( indicatorNum=5,
tacsatp,
minThreshold=0, # percentage of points to be included in the Minimal Convex Polygon
cellresX=0.05,
cellresY=0.025,
calcAreaMethod="Trapezoid",
plotMapTF=T,
exportTableName="",
exportGridName="")

cbind(ind5,ind5_2)
```

That was DCF indicator number 5, and changing the resolution made quite a change. Let's move on to indicator number 6 and 7.

```

library(adehabitat); library(gpclib)
#If it doesn't work:
# install.packages("adehabitat",repos=getOption("repos"))
# install.packages("gpclib",repos=getOption("repos"))

ind6 <- indicators(indicatorNum=6,
tacsatp,
pctThreshold=90,
cellresX=0.01,
cellresY=0.05,
calcAreaMethod="Trapezoid",
plotMapTF=TRUE,
exportTableName="",
exportGridName="")

#And now indicator number 7
ind7 <- indicators(indicatorNum=7,
tacsatp,
ltGear=c("TBB","OTB"),
inShapeArea="",
cellresX=0.01,
cellresY=0.05,
calcAreaMethod="Trapezoid",
minThreshold=0,
plotMapTF=TRUE,
exportGridName="")
```

For most of you these calculations will be taking place in a 'black box'. This doesn't mean you cannot use them.  Read the help files carefully and we encourange you to look at the code, which is the real beauty of R and Open Source.  If you want to compare indicator values across nations, where the assumptions on how to calculate aren't clear, this might serve you well enough.

However, if you want or need full control, you might be better off calculating these indicators by yourself. Find below some instructions on how to get started!

## Calculating indicators tailor-made ##
For a start, let's take a simple indicator to calculate ourselves. Let's define: the total surface area trawled by year. What we need to do for this exercise is:
  1. clean the tacsat and eflalo data
  1. merge them together
  1. perform activity analyses and take fishing activity only
  1. select bottom trawling gear only
  1. define the width of the gear
  1. define the grid

We have the choice to interpolate the data or not, but to keep things simple, we will not do that here. Obviously, the function to add a width to a gear, seen before, could be useful. However, we take less of a true spatial approach here (for most spatial analyses, your data structure is fine as vectors rather than as a grid!), whereby the 'addWidth' function would simply take too long to run.

```

data(tacsat);data(eflalo)

#Now clean the tacsat and eflalo dataset...
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]

#Again, the quick and dirty fishing activity analyses
tacsatp <- filterTacsat(tacsatp,st=c(2,6),hd=NULL,remDup=T)

#Select gears
table(tacsatp$LE_GEAR)
gears <- c("OTB","PTB","TBB")
tacsatp <- subset(tacsatp,LE_GEAR%in% gears)
gearWidth <- data.frame(LE_GEAR=gears,LE_WIDTH=c(87,2*87,24))

#Define the grid & subset the area
spatBound <- list(x=c(-5,10),y=c(48,62))
idxlon    <- which(tacsatp$SI_LONG >= spatBound$x[1] & tacsatp$SI_LONG <= spatBound$x[2])
idxlat    <- which(tacsatp$SI_LATI >= spatBound$y[1] & tacsatp$SI_LATI <= spatBound$y[2])
tacsatp   <- tacsatp[idxlon[which(idxlon %in% idxlat)],]
grd <- createGrid(xrange=c(unlist(spatBound$x)),
yrange=c(unlist(spatBound$y)),
resx=0.1,resy=0.05) #Let's take a tenth of an ICES rectangle for a start
```

So far so good, you've all seen these preparation steps a number of times now, but it's good to learn them by heart as you will need to go over these time and time again (look at [Practicals2](Practicals2.md)).

The next step is to start calculating the area trawled for each VMS ping. It might be handy here to store the results by year in some sort of data-frame again, so you can access it later, and also derive other results from it.

Area trawled = duration of trawling x width of gear x speed of a vessel.
Let's work this equation out for each of the VMS pings in the dataset.

```

tacsatp <- sortTacsat(tacsatp)

#Calculate duration of trawling per ping
tacsatp <- intervalTacsat(tacsatp,level="trip",fill.na=TRUE)
tacsatp$SI_INTV <- tacsatp$INTV
gearWidth$LE_GEAR <- ac(gearWidth$LE_GEAR)
tacsatp <- merge(tacsatp,gearWidth,by="LE_GEAR")
tacsatp$TR_AREA <- (tacsatp$SI_INTV / 60) * (tacsatp$LE_WIDTH / 1000) * (tacsatp$SI_SP *1.852)
```

For the full tacsat dataset we have now calculated the area trawled per ping, however, not yet per year. We could easily split up the tacsat dataset into the two years (1800 and 1801) and aggregate over the trawled area column. But, as indicated before, if we want a bit more flexibility, it might be good to store the data first and perform calculations later-on.

```

#Turn the defined grid into a spatial dataframe
sg <- SpatialGrid(grid=grd)
spgDF <- as(sg,"SpatialGridDataFrame")
slotNames(spgDF)

#We can simply add a dataframe to this object (but make sure it has the same size as all coordinates considered!
spgDF@data <- data.frame("Y1800"=rep(0,nrow(coordinates(spgDF))),
"Y1801"=rep(0,nrow(coordinates(spgDF))),
surface=rep(0,nrow(coordinates(spgDF))))
head(spgDF@data)

#Now align the VMS pings to the grid
#We have to turn the coordinates of the tacsat dataset into spatial points first
for(iYr in c(1800:1801)){
subTacsat <- subset(tacsatp,year(SI_DATIM)==iYr)
sp <- SpatialPoints(coordinates(subTacsat[,c("SI_LONG","SI_LATI")]))

#The spatial grid dataframe is a long list (nrow(spgDF@data)) and the idx tells me which row number (associated with a grid cell) the VMS ping is located in.
idx <- overlay(sg,sp)
spgDF@data[ac(aggregate(subTacsat$TR_AREA,          by=list(idx),FUN=sum,na.rm=T)[,1]), paste("Y",iYr,sep="")]    <- aggregate(subTacsat$TR_AREA,  by=list(idx),FUN=sum,na.rm=T)[,2]
}

#Have a look at some data (where Y1800 is not zero)
head(spgDF@data[which(spgDF@data$Y1800>0),])

#And calculate the indicator value
colSums(spgDF@data[,c("Y1800","Y1801")])
```

However, we just set this system up to be a bit more flexible and use the data for other questions too. Can we now easily move from total area trawled to percentage of the area untrawled? Or could we even see which area has been trawled in 1800 but not in 1801?

```

#Add surface of each gridcell to the dataframe
spgDF@data$surface <- surface(spgDF)$cellArea

#Percentage trawled
idxun <- which(spgDF@data[,"Y1800"] == 0)
sum(spgDF@data$surface[idxun]) / sum(spgDF@data$surface)

#Area trawled in 1800 but not in 1801
idx <- which(spgDF@data[,"Y1800"] > 0 & spgDF@data[,"Y1801"] == 0)
sum(spgDF@data$surface[idx]) / sum(spgDF@data$surface)
```

All in all, a lot to play with in a rather flexible way.

### Exercise 1 ###
  1. Calculate the total area trawled by the TBB gears only (in 1800)
  1. Calculate the proportion of the total surface the TBB has trawled (in 1800)
  1. In how many grid cells was the trawling intensity > 1? (i.e. where has the area trawled exceeded the surface of the area)