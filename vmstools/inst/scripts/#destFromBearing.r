library(vmstools)

data(tacsat)

  #Sort the VMS data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Filter the Tacsat data
tacsat     <- filterTacsat(tacsat,c(4,8),hd=NULL,remDup=T)

  #Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,res=100,method="cHs",params=list(fm=0.5,distscale=20,sigline=0.2,st=c(2,6)),headingAdjustment=0)

xrange        <- range(unlist(lapply(interpolation,function(x){return(range(x[-1,1]))})),na.rm=T)
yrange        <- range(unlist(lapply(interpolation,function(x){return(range(x[-1,2]))})),na.rm=T)
plot(interpolation[[1]][-1,1],interpolation[[1]][-1,2],type="l",pch=19,lwd=1,asp=1/lonLatRatio(interpolation[[1]][2,1],interpolation[[1]][2,2])[1],xlim=xrange,ylim=yrange,xlab="Longitude",ylab="Latitude")
for(i in 2:length(interpolation)){
  lines(interpolation[[i]][-1,1],interpolation[[i]][-1,2],type="l",pch=19,lwd=1,asp=1/lonLatRatio(interpolation[[i]][2,1],interpolation[[i]][2,2])[1])
}

interpolationGearWidth <- addWidth(interpolation,gearWidth=0.5)
plot(interpolationGearWidth,border="grey",col="grey",asp=1/lonLatRatio(interpolation[[1]][2,1],interpolation[[1]][2,2])[1],xlim=xrange,ylim=yrange,xlab="Longitude",ylab="Latitude")
box(); axis(1); axis(2); mtext(side=1,"Longitude",line=3); mtext(side=2,"Latitude",line=3)

  #Plot the interpolation
int <- 112
plot(interpolation[[int]][-1,1],interpolation[[int]][-1,2],type="l",pch=19,lwd=2,asp=1/lonLatRatio(interpolation[[int]][2,1],interpolation[[int]][2,2])[1])
xs <- interpolation[[int]][-1,1]
ys <- interpolation[[int]][-1,2]

  #Calculate the bearing towards and away from each point
bear1 <- bearing(xs[1:(length(xs)-2)],ys[1:(length(xs)-2)],xs[2:(length(xs)-1)],ys[2:(length(xs)-1)])
bear2 <- bearing(xs[2:(length(xs)-1)],ys[2:(length(xs)-1)],xs[3:length(xs)],ys[3:length(xs)])
avbear<- atan2(mapply(sum,sin(bear1*(pi/180))+sin(bear2*(pi/180))),mapply(sum,cos(bear1*(pi/180))+cos(bear2*(pi/180))))*(180/pi)

  #Take the average of the two
avbear<- c(avbear[1],avbear,avbear[length(avbear)])

  #Calculate the destinated point taking a begin point, a bearing and a certain distance to travel
outpointr <- destFromBearing(xs,ys,(avbear+90+360)%%360,0.5)
outpointl <- destFromBearing(xs,ys,(avbear-90+360)%%360,0.5)

  #Plot these lines
lines(outpointr[,1],outpointr[,2],col="red",lty=2)
lines(outpointl[,1],outpointl[,2],col="red",lty=2)


  #Create polygons from it
for(i in 1:(nrow(outpointr)-1)){
  polygon(x=c(outpointr[i,1],outpointl[i,1],outpointl[i+1,1],outpointr[i+1,1]),y=c(outpointr[i,2],outpointl[i,2],outpointl[i+1,2],outpointr[i+1,2]),col="black")
}

pols <- list()
for(i in 1:(nrow(outpointr)-1)){
  pols[[i]] <- Polygon(cbind(c(outpointr[i,1],outpointl[i,1],outpointl[i+1,1],outpointr[i+1,1],outpointr[i,1]),
                             c(outpointr[i,2],outpointl[i,2],outpointl[i+1,2],outpointr[i+1,2],outpointr[i,2])))
}

polys <- list()
for(i in 1:10){
  polys[[i]] <- Polygons(pols,ID=ac(i))
}
spPolys <- SpatialPolygons(polys)
plot(spPolys,col="black")





  #Destination:
  #This function takes a starting x,y position, a bearing and a distance to follow the initial bearing.
  #To get this bearing, you have to compute the bearing from an towards the x,y position. +/- 90 degrees
  #as this the most outer point taken from the point you are heading in. Then travel along that heading
  #for a certain km's and you'll get to your endpoint.
  
