library(vmstools)

data(tacsat)

  #Sort the VMS data
tacsat     <- sortTacsat(tacsat)
tacsat     <- tacsat[1:1000,]

  #Filter the Tacsat data
tacsat     <- filterTacsat(tacsat,c(4,8),hd=NULL,remDup=T)

  #Interpolate the VMS data
interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,res=100,method="cHs",params=list(fm=0.5,distscale=20,sigline=0.2,st=c(2,6)),headingAdjustment=0)

  #Plot the interpolation
int <- 120
plot(interpolation[[int]][-1,1],interpolation[[int]][-1,2],type="l",pch=19,lwd=2,asp=1/lonLatRatio(interpolation[[int]][2,1],interpolation[[int]][2,2])[1])
xs <- interpolation[[int]][-1,1]
ys <- interpolation[[int]][-1,2]

  #Calculate the bearing towards and away from each point
bear1 <- bearing(xs[1:(length(xs)-2)],ys[1:(length(xs)-2)],xs[2:(length(xs)-1)],ys[2:(length(xs)-1)])
bear2 <- bearing(xs[2:(length(xs)-1)],ys[2:(length(xs)-1)],xs[3:length(xs)],ys[3:length(xs)])
avbear<- (bear1+bear2)/2
  #Take the average of the two
avbear<- c(avbear[1],avbear,avbear[length(avbear)])

  #Calculate the destinated point taking a begin point, a bearing and a certain distance to travel
outpointr <- destFromBearing(xs,ys,(avbear+90)%%360,0.5)
outpointl <- destFromBearing(xs,ys,(avbear-90)%%360,0.5)

  #Plot these lines
lines(outpointr[,1],outpointr[,2],col="red",lty=2)
lines(outpointl[,1],outpointl[,2],col="red",lty=2)



plot(interpolation[[int]][-1,1],interpolation[[int]][-1,2],type="l",pch=19,lwd=2,asp=1/lonLatRatio(interpolation[[int]][2,1],interpolation[[int]][2,2])[1])

rxs <- range(xs); rys <- range(ys)
points(rxs,rys,pch=19,col="darkgreen")
points(rev(rxs),rys,pch=19,col="darkgreen")
lines(rxs,rys,col="green")
lines(rev(rxs),rys,col="green")

  #Destination:
  #This function takes a starting x,y position, a bearing and a distance to follow the initial bearing.
  #To get this bearing, you have to compute the bearing from an towards the x,y position. +/- 90 degrees
  #as this the most outer point taken from the point you are heading in. Then travel along that heading
  #for a certain km's and you'll get to your endpoint.
  
destFromBearing <- function(lon1,lat1,bearing,distance){

                  pd    <- pi/180
                  dist  <- (distance/6371)
                  y1    <- lat1 * pd
                  x1    <- lon1 * pd
                  bear  <- bearing * pd
                  
                  y2  <- asin(sin(y1) * cos(dist) + cos(y1) * sin(dist) * cos(bear))
                  x2  <- x1 + atan2(sin(bear) * sin(dist) * cos(y1),cos(dist) - sin(y1) * sin(y2))

                  x2  <- (x2 + 3*pi) %% (2*pi) - pi
                  y2  <- y2 / pd
                  x2  <- x2 / pd
                return(cbind(x2,y2))}