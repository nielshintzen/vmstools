plotSingleInterpolation <- function(interpolation,reference){
                              intIdx  <- interpolation[1,]
                              xrange  <- range(c(interpolation[-1,1],reference[intIdx[1]:intIdx[2],"declon"]),na.rm=T)
                              yrange  <- range(c(interpolation[-1,2],reference[intIdx[1]:intIdx[2],"declat"]),na.rm=T)

                                          plot(interpolation[-1,1],interpolation[-1,2],xlim=xrange,ylim=yrange,type="l",
                                               xlab="Longitude",ylab="Latitude",lwd=2)
                                          lines(reference[intIdx[1]:intIdx[2],"declon"],reference[intIdx[1]:intIdx[2],"declat"],
                                               col="red",lwd=2,lty=2)
                                          }