plotInterpolation(interpolation,buffer=T,bufferSize=1)

plotInterpolation <- function(interpolation,buffer=F,bufferSize=0){

                          xrange <- range(unlist(lapply(interpolation,function(x){return(x[-1,1])})),na.rm=T)
                          yrange <- range(unlist(lapply(interpolation,function(x){return(x[-1,2])})),na.rm=T)

                          plot(interpolation[[1]][-1,1],interpolation[[1]][-1,2],type="l",lwd=2,xlim=range(pretty(xrange)),ylim=range(pretty(yrange)),asp=lonLatRatio((xrange[2]-xrange[1])/2+xrange,(yrange[2]-yrange[1])/2+yrange)[1])
                          if(length(interpolation) > 1){
                            for(int in 2:length(interpolation)){

                              xs <- interpolation[[int]][-1,1]
                              ys <- interpolation[[int]][-1,2]
                              lines(xs,ys,lwd=2)
                              if(buffer == T & bufferSize > 0){
                                  #Calculate the bearing towards and away from each point
                                bear1 <- bearing(xs[1:(length(xs)-2)],ys[1:(length(xs)-2)],xs[2:(length(xs)-1)],ys[2:(length(xs)-1)])
                                bear2 <- bearing(xs[2:(length(xs)-1)],ys[2:(length(xs)-1)],xs[3:length(xs)],ys[3:length(xs)])
                                avbear<- (bear1+bear2)/2
                                  #Take the average of the two
                                avbear<- c(avbear[1],avbear,avbear[length(avbear)])

                                  #Calculate the destinated point taking a begin point, a bearing and a certain distance to travel
                                outpointr <- destFromBearing(xs,ys,(avbear+90)%%360,bufferSize)
                                outpointl <- destFromBearing(xs,ys,(avbear-90)%%360,bufferSize)

                                  #Plot these lines
                                lines(outpointr[,1],outpointr[,2],col="red",lty=2)
                                lines(outpointl[,1],outpointl[,2],col="red",lty=2)
                              }
                            }
                          }
                    }




