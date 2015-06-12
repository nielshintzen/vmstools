# Practical 5: Interpolation methods for VMS tracks #

## Introduction ##
Sometimes, even scaling up from the 'low' spatial resolution of logbooks (ICES rectangles) to the more detailed VMS resolution is not enough.  In this case you may want to interpolate the fishing tracks between your VMS data points, usually separated by 2 hours, and create a number of 'intermediate' points.  We call these 'interpolated tracks'.

Before you start to consider this type of analysis, you have to be aware of a number of issues.
  1. You are now 'making up' data points to a certain extent
  1. You might associate these interpolated points with catches, so keep an eye on reality
  1. Most often you only want to interpolate between fishing pings, so your analyses of activity is very important!
  1. Is your computer up for the test, because you will expand your dataset even more (and this will require even more memory and storage)

That should be enough warnings to convince you how useful these things actually are, and please, do try this at home, but be aware of the assumptions when you present something to, e.g. stakeholders or a government ministry.

## Interpolate VMS data ##
In VMStools we only distinguish between two types of interpolation:
  1. Linear interpolation
  1. Cubic Hermite Spline interpolation (according to Hintzen et al. 2010, Fisheries Research)

Although the technicalities might be a bit too much for now, the practical interpretation is actually quite straightforward.

```

data(tacsat)

#Sort the data and make a subset
tacsat <- sortTacsat(tacsat)[1:1000,]
tacsat <- filterTacsat(tacsat,st=c(2,6),hd=NULL,remDup=T)

#Do not be put off by the large number of arguments here
#Interpolate according to a straight line interpolation
interpolationSL <- interpolateTacsat(tacsat,interval=120,margin=10,res=100,method="SL",params=list(fm=0.5,distscale=20,sigline=0.2,st=c(2,6)),headingAdjustment=0)

#Interpolate according to the cubic-hermite spline interpolation
interpolationcHs <- interpolateTacsat(tacsat,interval=120,margin=10,res=100,method="cHs",params=list(fm=0.2,distscale=20,sigline=0.2,st=c(2,6)),headingAdjustment=0)

#The returned interpolation is an R 'list', let's look at the first interpolation
plot(interpolationSL[[1]][-1,],type="l",xlim=c(3.17,3.24))
points(interpolationSL[[1]][c(2,101),],pch=19,col="blue")
lines(interpolationcHs[[1]][-1,],col="red")
```

Ok, we've managed to interpolate some of our data but interpolations were not possible in every instances, as we only got:

```

#How many interpolations do I have?
length(interpolationSL)
#How many tacsat entries do I have?
nrow(tacsat)
```

216 interpolations, while we had 414 tacsat points. In theory (if they are all of the same vessel) we could have had 413 interpolations.

The format of the interpolation allows some 'lapply' functions to be executed directly on the object.  In other cases, however, it might be easier to convert the interpolated data back to the default tacsat format.

```

#Convert the interpolation to tacsat style data
tacsatInt <- interpolation2Tacsat(interpolationcHs,tacsat,npoints=10)

#Each interpolation is split into 10 equally spaced points (including the start and the end point). In total 8 new points are added per interpolation.
nrow(tacsatInt)
nrow(tacsat)
(nrow(tacsatInt)-nrow(tacsat))/8 #I've seen that number before

#What is the distance traveled between each of these interpolations?
sum(distanceInterpolation(interpolationSL))
sum(distanceInterpolation(interpolationcHs))
```

The only way to really test if either of these two interpolations is any better is to test it against a much higher resolution dataset.

```

#Make use of a 6-minute interval dataset, and subset that every 120 minutes
#Load some high frequency VMS data
data(VMShf)

#And go through some formatting issues and very basic cleaning...
colnames(VMShf) <- c("VE_REF","SI_LATI","SI_LONG","SI_SP","SI_HE","SI_DATIM")
VMShf$VE_COU <- "Atlantis"
VMShf$SI_DATE<- format(VMShf$SI_DATIM,format="%d/%m/%Y")
VMShf$SI_TIME<- format(VMShf$SI_DATIM,format="%H:%M")
VMShf <- formatTacsat(VMShf)
VMShf <- sortTacsat(VMShf)
VMShf <- filterTacsat(VMShf,st=c(2,6),hd=NULL,remDup=T)

#cHs first
cHs <-interpolateTacsat(VMShf,interval=120,margin=10,res=100,method="cHs",params=list(fm=0.2,distscale=20,sigline=0.2,st=c(2,6)),headingAdjustment=1)

#Then straight line
SL <-interpolateTacsat(VMShf,interval=120,margin=10,res=100,method="SL",params=list(fm=0.2,distscale=20,sigline=0.2,st=c(2,6)),headingAdjustment=1)

#Visualize the results, note the diffInter function, what does it do?
plot(diffInter(cHs,VMShf)[,"mean"],pch=19,ylab="difference",ylim=c(0,10))
abline(h=1,col=2)
points(diffInter(SL,VMShf)[,"mean"],col="blue",pch=19)
length(which(diffInter(cHs,VMShf)[,"mean"] < diffInter(SL,VMShf)[,"mean"]))
sum(diffInter(cHs,VMShf)[,"mean"])
sum(diffInter(SL,VMShf)[,"mean"])
```

Obviously, we could do a bit of tuning here and there to make the performance of the interpolation even better, as this was just a quick and dirty example showing the benefits.

### Exercise 1 ###
  1. Plot the third interpolation, why is the curve so different from the straight line?
  1. What is the difference in distance traveled between the two types of interpolation? (try and use the 'distance' function to do it yourself)
  1. What would be the associated speed of the vessel under the two different interpolation scenarios?

## From here onwards ##

Once the tacsat data is interpolated some other possibilities remain. It is, for example, possible to add a certain width of a gear to the interpolation. This might be handy when you want to reconstruct trawling tracks at approximately the right scale. Or you could focus on the uncertainty associated with interpolated tracks to test if areas have remained unfished for some time.

### Adding gear width to interpolation ###
Especially when your study area is very small (e.g. only a few kilometers), it is worth to start considering the width of the gear in fishing track pictures. One feature tested before was Google Earth pictures where, by zooming in, the actual with of the gear could be visualized. As we now move into rather sophisticated tools, I do show some of the underlying tools that are used, but wrapped, inside the functions.

First of all, select a vessel for which we know the gear width (e.g. a beam trawler that fishes two 12 meter beams).
```

data(eflalo); data(tacsat)
#Select one vessel
unique(subset(eflalo,LE_GEAR == "TBB"))$VE_REF[1]

subTacsat <- subset(tacsat,VE_REF == "238")

#Interpolate the dataset
cHs <-interpolateTacsat(subTacsat,interval=120,margin=10,res=100,method="cHs",params=list(fm=0.2,distscale=20,sigline=0.2,st=c(2,6)),headingAdjustment=0)

#Make a picture of all interpolations first
#Get the ranges of the total picture
ranges <- do.call(rbind,lapply(cHs,function(x){return(apply(x[-1,],2,range))}))
xrange <- range(ranges[,"x"])
yrange <- range(ranges[,"y"])

#Plot the results
plot(subTacsat$SI_LONG,subTacsat$SI_LATI,xlim=xrange,ylim=yrange,pch=19,col="blue",xlab="Longitude",ylab="Latitude")
for(iInt in 1:length(cHs)){
lines(cHs[[iInt]][-1,1],cHs[[iInt]][-1,2])}
```

What we now basically would like to do is something close to this:

```

#Change width of the line to 2 instead of 1
for(iInt in 1:length(cHs)){
lines(cHs[[iInt]][-1,1],cHs[[iInt]][-1,2],lwd=2)}
```

That did give it more wide gear tracks, but how wide are they really? The total distance from the left hand side of the picture to the right hand side equals approximately:

```

distance(3.6,51.8,4.5,51.8) #~62km
```

Nearly 62 km! So, using a thicker line to plot is not really accurate. Let's do it the 'right' way.

```

interWidth <- addWidth(cHs,gearWidth=0.024)
x11(); plot(interWidth); map.axes()
points(subTacsat$SI_LONG,subTacsat$SI_LATI,col="blue",pch=19)
```

Again, you can question how useful this really is if you look at a picture with dimension 60x70km... It would be better to zoom in a bit more and have another look and compare the results with a 'normal' line drawing of R.

```

#Let's look at 2km width, how many degrees is that?
distance(3.6,51.8,3.7,51.8) #0.1 degree longitude = 6.8km, so 1 km is 0.015 degree longitude.
distance(3.6,51.8,3.615,51.8)

#In latitude, approx half of that step size
distance(3.6,51.8,3.6,51.81)

#Take a VMS point and take 100 meter on each side of it
plot(interWidth,ylim=c(51.83098 - (0.0075),51.83098 + (0.0075 )),xlim=c(3.778870 - (0.015),3.778870 + (0.015)))
points(subTacsat$SI_LONG,subTacsat$SI_LATI,pch=19)
map.axes()
for(iInt in 1:length(cHs)){
lines(cHs[[iInt]][-1,1],cHs[[iInt]][-1,2],col="red")} #The interpolated line only
```

What you can see in the picture is the VMS point (the dot), the interpolated line (without gear width, in red) and the interpolated line with gear width in black squares. What happens under the hood is that for each little step taken in the interpolation (each interpolation consists of e.g. 100 small steps i.e. new points), it's bearing is calculated towards the next small step. Having this knowledge, the width of the gear must extend in a perpendicular direction (90 and -90 degrees from the heading of the vessel). Keeping this bearing for a distance equal to the width of the gear gives us the outer points of the black boxes. All black boxes are thereafter combined to give a fishing track. So, the function bearing plays an important role here, you might want to use it for other purposes too.

```

bearing(3.7,51.8,3.6,51.8) #heading into a 270 degree compass course

#If you are interested where you will end up with a certain bearing and distance to travel, you can use this
distance(3.7,51.8,3.6,51.8) #what distance have I traveled?
destFromBearing(3.7,51.8,270.0393,6.876387) #starting from the same point, traveling with the bearing and distance calculated, do I end up at the same point?
```

It is recommended to be specific when calculating the area trawled using the 'addWidth' functionality as it requires a very large number of calculations which will most certainly slow down your analyses.

### Calculating uncertainty of interpolation ###
'Uncertainty' is an idea common to most scientists. This uncertainty can, however, also provide interesting results, and be used in, e.g. indicator calculations!

Obviously, interpolation isn't without uncertainty either. Hence, tools have been developed to show / calculate this uncertainty and to allow us to do other analyses on the basis of those results. It should be stressed, however, that we are talking about uncertainty in the track of the interpolation, not the VMS position itself.  Remember that the VMS 'pings' aren't without error (although it is probably negligible in the grand scheme of things).

From the interpolation and VMS points themselves, we can already distill a bit of guidance on how to construct 'confidence intervals'.
  1. Confidence at the VMS pings equals 1
  1. The further away from, either of these two points, the lower the certainty.
  1. The further away from the interpolation, the lower the certainty (as we assume that the interpolation is the best possible track between the two pings)

We like to do these calculations on a grid again, as that allows us to plot the results too.
```

data(tacsat)

#Sort the Tacsat data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

#Filter the Tacsat data
tacsat          <- filterTacsat(tacsat,c(2,6),hd=NULL,remDup=T)

#Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
res=100,method="cHs",params=list(fm=0.5,distscale=20,
sigline=0.2,st=c(2,6)),headingAdjustment=0)

#Create the final grid where all interpolations should fit on
xrange        <- c(2.1,2.6); yrange <- c(51.3,51.65)
grid          <- createGrid(xrange,yrange,resx=0.001,resy=0.0005)

#Do some fancy stuff and calculate the Confidence Interval
res           <- calculateCI(interpolation[[4]],
tacsat[interpolation[[4]][1,],],
params=list(fm=0.25,distscale=3.1,sigline=0.4,st=c(2,6)),
grid=grid,
plot=TRUE)

#Scale all values between 0 and 1
res <- (res[[1]] - min(res[[1]])) / (max(res[[1]]) - min(res[[1]]))

levelplot(matrix(res,ncol=458,nrow=380)[,380:1])
```

That was the code to just calculate one confidence interval. But what if I want the confidence interval of a whole trip. That requires a short loop! Let's see how far we get.

```

data(tacsat)

#Sort the Tacsat data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

#Filter the Tacsat data
tacsat          <- filterTacsat(tacsat,c(2,6),hd=NULL,remDup=T)
#Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
res=100,method="cHs",params=list(fm=0.5,distscale=20,
sigline=0.2,st=c(2,6)),headingAdjustment=0)
```

Ok, data has been prepared and interpolated, now define the grid

```

ranges        <- do.call(rbind,lapply(interpolation,function(x){return(apply(x[-1,],2,range))}))
xrange        <- range(ranges[,1]); yrange <- range(ranges[,2])
grid          <- createGrid(xrange,yrange,resx=0.01,resy=0.005)
sPDF          <- createGrid(xrange,yrange,resx=0.01,resy=0.005,type="SpatialGridDataFrame")

#Now I don't add zeros, but ones instead, because it's about the chance of NOT being trawled which declines with the values of the CI
sPDF@data     <- data.frame(rep(1,nrow(sPDF@data)))
colnames(sPDF@data) <- "data"
```

I've created a data frame where I can store the chance that an area is untrawled, which obviously ranges between 0 and 1. I call this the sPDF (spatial data frame). Hereafter I need to fill the spatial dataframe based on the calculations of the CI. The functions returns the index of the dataframe where data on the CI should be stored (not in every grid cell, but just a selection) and the value of the CI cells.

```

library(maps); library(mapdata); library(RColorBrewer)
#Start the loop
for(iNt in 1:length(interpolation)){
print(iNt)
resCI <- calculateCI(interpolation[[iNt]],tacsat[interpolation[[iNt]][1,],],
params=list(fm=0.5,distscale=3.1,sigline=0.42,st=c(2,6)),grid=grid,
plot=F)

idx       <- getGridIndex(coordinates(resCI),grid,all.inside=F)
#some coordinates might be located outside of the grid, we need to check and correct for that
idxdat    <- which(is.na(idx)==F)
sPDF@data$data[idx[idxdat]] <- (1-resCI@data$data[idxdat]) * sPDF@data$data[idx[idxdat]]
}#end for loop
```

As I am only interested in those areas where I do have trawling, I simply only use that part where the chance of no trawling is smaller than 1. Thereafter I associate the value of the data with a color, which makes it easier to plot afterwards.

```

#Select those areas where trawling intensity is smaller than 1
idxwithdata <- which(sPDF@data<1)

#Get the coordinates of these areas
coords      <- data.frame(coordinates(sPDF)[idxwithdata,])
coords$data <- sPDF@data[idxwithdata,]

#some nasty code here just to determine in which category (in total 9) the data belongs. This category automatically becomes the number of the color to use
coords$color<- apply(abs(outer(seq(0,1,length.out=9),coords$data,"-")),2,which.min)

#Define the ranges in longitude and lattitude
xrange <- range(coords[,1]); yrange <- range(coords[,2])
map("worldHires",fill=T,col="darkgreen",xlim=xrange,ylim=c(yrange[1],54.6))
map.axes()
color <- rev(brewer.pal(9,"YlOrRd"))
```

Time to plot, a 'for' loop with a 'polygon' call is not the best way to do it, but it works for now.

```

for(i in 1:nrow(coords)){
polygon(x=c(coords[i,1],coords[i,1]+0.01,coords[i,1]+0.01,coords[i,1]),
y=c(coords[i,2],coords[i,2],
coords[i,2]+0.005,coords[i,2]+0.005),
col=color[coords[i,"color"]],lwd=1,border=NA)
}
#Some plotting occurred on top of the map, so let's add the map again (but do not overwrite, hence the 'add=T' argument
map("worldHires",fill=T,col="darkgreen",xlim=xrange,ylim=yrange,add=T)
```

That gives us a nice picture, let's see where the original VMS datapoints were located!

```

#compare to VMS points
x11()
map("worldHires",fill=T,col="darkgreen",xlim=xrange,ylim=c(yrange[1],54.6))
map.axes()
points(tacsat$SI_LONG,tacsat$SI_LATI,col=1,pch=19,cex=0.5)
```

That is pretty much all VMStools can help you with related to interpolations. Go ahead and have fun!