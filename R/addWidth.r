addWidth <- function(interpolation,gearWidth){


allPolygons <- list()
counter     <- 0
for(i in 1:length(interpolation)){

    #Take the interpolated points
  xs <- interpolation[[i]][-1,1]
  ys <- interpolation[[i]][-1,2]

  if(all(is.na(xs))==F & all(is.na(ys))==FALSE){
    counter <- counter + 1
      #Calculate the bearing towards and away from each point
    bear1 <- bearing(xs[1:(length(xs)-2)],ys[1:(length(xs)-2)],xs[2:(length(xs)-1)],ys[2:(length(xs)-1)])
    bear2 <- bearing(xs[2:(length(xs)-1)],ys[2:(length(xs)-1)],xs[3:length(xs)],ys[3:length(xs)])
    avbear<- atan2(mapply(sum,sin(bear1*(pi/180))+sin(bear2*(pi/180))),mapply(sum,cos(bear1*(pi/180))+cos(bear2*(pi/180))))*(180/pi)

      #Take the average of the two
    avbear<- c(avbear[1],avbear,avbear[length(avbear)])

      #Calculate the destinated point taking a begin point, a bearing and a certain distance to travel
    outpointr <- destFromBearing(xs,ys,(avbear+90+360)%%360,gearWidth/2)
    outpointl <- destFromBearing(xs,ys,(avbear-90+360)%%360,gearWidth/2)

    singlePolygons <- list()
    for(j in 1:(nrow(outpointr)-1)){
      singlePolygons[[j]]   <- Polygon(cbind(c(outpointr[j,1],outpointl[j,1],outpointl[j+1,1],outpointr[j+1,1],outpointr[j,1]),
                                             c(outpointr[j,2],outpointl[j,2],outpointl[j+1,2],outpointr[j+1,2],outpointr[j,2])))
    }
    allPolygons[[counter]] <- Polygons(singlePolygons,ID=ac(counter))
  }
}
return(SpatialPolygons(allPolygons))}
