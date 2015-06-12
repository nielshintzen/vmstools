# Practical 4: How to link VMS and logbook data to 'spatial' grids #

## Introduction ##
In general, logbook data are recorded at the ICES rectangle level. VMS data, however, are not associated with any spatial scale (other than taking account of the uncertainty in the data themselves which may be a limiting factor). In many instances, however, the output of any analysis must be delivered at some spatial scale. If that spatial scale is ICES rectangle for logbook data, you're lucky, but most often it requires a bit more work.

Here we will provide some tools to help you out with such jobs.

## Converting between spatial notations ##
What if you are required to report VMS data analyses on the ICES rectangle? Will you go through some complicated procedure to find the closest ICES rectangle to your data, or just use the embedded functionality in VMStools? Let's see how that would work.

```

#Start afresh!
library(vmstools)
rm(list=ls()) #This removes all objects from R

#load the data
data(tacsat); data(eflalo)

#Add ICES rectangles to your tacsat data, based on longitude and latitude.
tacsat$LE_RECT <- ICESrectangle(tacsat)

#That was easy, can we also convert back?
pos <- ICESrectangle2LonLat(tacsat$LE_RECT,midpoint=T)

tacsatCompare <- cbind(tacsat,pos)

head(tacsatCompare)
>    VE_COU VE_REF  SI_LATI  SI_LONG    SI_DATE  SI_TIME SI_SP SI_HE LE_RECT SI_LATI SI_LONG
658237 Atlantis     10 51.44772 3.583731 04/05/1800 13:32:00     0    25    31F3   51.25     3.5
658238 Atlantis     10 51.44067 3.583847 04/05/1800 15:28:00     0    25    31F3   51.25     3.5
662888 Atlantis     10 51.44074 3.595529 04/05/1800 17:22:00     0    25    31F3   51.25     3.5
662889 Atlantis     10 51.44315 3.586862 04/05/1800 19:18:00     0    25    31F3   51.25     3.5
662890 Atlantis     10 51.44351 3.586948 04/05/1800 21:12:00     0    25    31F3   51.25     3.5
662891 Atlantis     10 51.44811 3.593958 04/05/1800 23:08:00     0    25    31F3   51.25     3.5
```

### Exercise 1 ###
  1. When executing the function 'ICESrectangle2LonLat' a warning was given. What does this warning mean?
  1. Can you find out the locations of the points not converted back to longitude,latitude? Does this make sense?

We can of course do similar things with the logbook data

```

pos <- ICESrectangle2LonLat(eflalo$LE_RECT,midpoint=T)
eflalop <- cbind(eflalo,pos)

#Let's find the ICES area associated with the ICES rectangle
#We need the longitude, latitude position though from the ICES rectangle midpoint from the step before
data(ICESareas)
eflalop$LE_AREA <- ICESarea(eflalop,ICESareas)
```

Obviously, we can do the same for the VMS data, so try it yourself! One other type of conversion available is the spatial <b>C-square</b> notation. This notation was developed in Australia, and has many handy features, ie. it can be used for the entire globe.  It is now implemented as the default FishFrame spatial scale. The C-square notation is flexible since, although it is a standard 'naming convention', it can also take many different resolutions. Let's convert some VMS data into C-square notation.

```

#For this example, take the first six positions only
somePos <- head(tacsat[,c("SI_LONG","SI_LATI")])
CSquare(somePos$SI_LONG,somePos$SI_LATI,degrees=10)
CSquare(somePos$SI_LONG,somePos$SI_LATI,degrees=10)
CSquare(somePos$SI_LONG,somePos$SI_LATI,degrees=5)
CSquare(somePos$SI_LONG,somePos$SI_LATI,degrees=0.5)
CSquare(somePos$SI_LONG,somePos$SI_LATI,degrees=0.1)
CSquare(somePos$SI_LONG,somePos$SI_LATI,degrees=0.05)
CSquare(somePos$SI_LONG,somePos$SI_LATI,degrees=0.01)
```

The length of the CSquare notation becomes longer and longer the higher the resolution gets.

## Analyse on a defined spatial scale ##

So far, we've only converted VMS and logbook data into different spatial notations. A next step would be to analyse the data based on the converted spatial scale.

```

#Let's calculate the catches of plaice in the total eflalo dataset at the ICES rectangle level
aggregate(eflalop$LE_KG_PLE,by=list(eflalop$LE_RECT),FUN=sum,na.rm=T)

#This would be the picture associated with that catch distribution
plotTools(eflalop,level="ICESrectangle",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,log=FALSE,color=NULL,control.tacsat=list(clm=NULL),control.eflalo=list(clm=c("LE_KG_PLE")))

#Now at a slightly finer scale, the VMS data
tacsatp <- mergeEflalo2Tacsat(eflalo,tacsat)
tacsatp$IDX       <- 1:nrow(tacsatp)

#A quick and dirty rule which says that all speeds between 1 and 6 are fishing
tacsatFilter      <- filterTacsat(tacsatp,st=c(1,6),hd=NULL,remDup=T)
tacsatp$SI_STATE  <- 0
tacsatp$SI_STATE[tacsatFilter$IDX] <- 1

#I needed the activity distinction because I only want to dispatch catches or values over fishing pings!
tacsatEflalo      <- splitAmongPings(tacsat=tacsatp,eflalo=eflalo,
variable="all",level="day",conserve=T)

#Let's count the number of rows in tacsatEflalo, and compare with the original tacsatp file
nrow(tacsatp)
nrow(tacsatEflalo)
table(tacsatp$SI_STATE) #does this make sense?

#Add a grid notation
tacsatEflalo$LE_SPAT <- CSquare(tacsatEflalo$SI_LONG,tacsatEflalo$SI_LATI,degrees=0.1)
aggregate(tacsatEflalo$LE_KG_PLE,by=list(tacsatEflalo$LE_SPAT),FUN=sum,na.rm=T)

#In a picture, it looks like this
x11()
plotTools(tacsatEflalo,level="gridcell",xlim=c(-5,10),ylim=c(48,62),zlim=NULL,log=F,gridcell=c(0.1,0.1),color=NULL,control.tacsat=list(clm="LE_KG_PLE"))
```

Observe the difference! This is a much more detailed picture of plaice catch distribution. Note the definition of the Plaice Box which was gradually closed to large beamers between 1989 and 1994.

### Exercise 2 ###
  1. Aggregate COD catches on a CSquare notation of 0.5 degrees
  1. If you are handy with R, can you make your own plot? (else, use plotTools to do it for you)

## Define your own grid ##
Although the spatial grids just discussed might be suitable in a large number of occasions, there are still some instances where you need to define your own grid. Within the core developers group of VMStools, we have decided not to use any km or nm based grids but only longitude and lattitude based grids. As all calculations are based on the assumption that the earth is spherical, we wouldn't be able to use these tools anymore if we would move to km or nm notation. We would be happy to embed successful code to perform these analyses in the package.

VMStools relies on a secondary package for spatial data storage and analyses, the sp package by Bivand et al. If you plan to do a lot of spatial analyses in R, in combination with VMS (basically the type of analyses you'll see from here and onwards) I can recommend you read the first chapters of his book <i><a href='http://books.google.nl/books?hl=nl&lr=&id=4gg8yx_lcj0C&oi=fnd&pg=PR7&dq=bivand,+sp,+R&ots=yzdH3v5-7G&sig=cGG9UC3rlreQxzlVROJXTKbxmyE#v=onepage&q=bivand%2C%20sp%2C%20R&f=false'>Applied spatial analyses with R (2008)</a></i>, especially chapter 2 on the different classes. It is not an easy read, but it will help you to understand how most of these analyses work together with polygons, shapefiles and data storage.

Let's define our first own grid

```

#Define your ranges in longitude (x) and latitude (y)
x <- seq(-5,5,length.out=10)
y <- seq(50,55,length.out=10)
xrange <- range(x,na.rm=T)
yrange <- range(y,na.rm=T)

#returns a grid with 11 cells in the x-direction and 11 in the y-direction
Grid <- createGrid(xrange,yrange,1,0.5,type="SpatialPixels") #just like an ICES rectangle grid

#Let's see what is stored in the 'Grid' object
slotNames(Grid)
>[1] "coords"      "grid"        "grid.index"  "bbox"        "proj4string"
Grid@coords #is the same as...
coordinates(Grid) #I prefer this way, works more often with other 'sp' objects too!

plot(Grid); map.axes()
```

Now you've defined your own grid, but obviously you want to link, e.g. your VMS data points, to this grid. In other words, you want to know in which grid cell each of your VMS positions would fall. Why is this handy? Primarily because you can now calculate certain indicators based on a user-defined grid (e.g. kg of cod in a certain area) and it allows you to combine different data sources within the same grid! But more on that later... First, let's see in which grid cell all theVMS positions fall.

```

#Get some tacsat positions, and link them to the Grid
somePos <- head(tacsat[,c("SI_LONG","SI_LATI")])
points(SpatialPoints(coordinates(somePos)),col="red",pch=19)
overlay(Grid,SpatialPoints(coordinates(somePos)))
>[1] 87 87 87 87 87 87

#This output really means that all 6 positions in 'somePos' are located in grid cell 87 of coordinates(Grid)
coordinates(Grid)[87,]

#Now for the full tacsat dataset, and summarize the result immediately
table(overlay(Grid,SpatialPoints(coordinates(tacsat[,c("SI_LONG","SI_LATI")]))))
```

### Exercise 3 ###
  1. Which grid cell (in coordinate notation) is associated with the highest number of VMS records?
  1. Plot it onto the grid as in the example above

All you actually need to aggregate again over a certain grid is the 'overlay' function, and the definition of your own grid. Once you get the results of the overlay back, you can use this to summarize your results. An example:

```

#Turn tacsat in the 'coordinates' class
coordTacsat <- coordinates(tacsat[,c("SI_LONG","SI_LATI")])

#And thereafter turn it into the SpatialPoints class
SPtacsat <- SpatialPoints(coordTacsat)

#Which grid cells do match with the tacsat points?
index <- overlay(Grid,SPtacsat) #the overlay function is being replaced by 'over', overlay still works though.
tacsat$INDEX <- index

#Do something whacky for a change, and take the mean heading in a grid cell
aggregate(tacsat$SI_HE,by=list(tacsat$INDEX),FUN=mean,na.rm=T)
```