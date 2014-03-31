#- Calculate proportion of land in ICESrectangles

resx    <- 1
resy    <- 0.5
xrange  <- seq(-4.5,10.5,resx)
yrange  <- seq(47.75,62.25,resy)

rects   <- expand.grid(SI_LONG=xrange,SI_LATI=yrange)
rects$LE_RECT <- ICESrectangle(rects)

#- Convert rects to SpatialPolygons
rectsSP <- lonLat2SpatialPolygons(lst=lapply(as.list(1:nrow(rects)),function(x){
              data.frame(SI_LONG=c(rects$SI_LONG[x]-0.5 ,rects$SI_LONG[x]+0.5,
                                   rects$SI_LONG[x]+0.5 ,rects$SI_LONG[x]-0.5),
                         SI_LATI=c(rects$SI_LATI[x]-0.25,rects$SI_LATI[x]-0.25,
                                   rects$SI_LATI[x]+0.25,rects$SI_LATI[x]+0.25))}))

#- Load Europa and turn it into SpatialPolygons too
data(europa)
eurSP   <-  lonLat2SpatialPolygons(lst=lapply(as.list(sort(unique(europa$SID))),
                 function(x){data.frame(SI_LONG=subset(europa,SID==x)$X,
                                      SI_LATI=subset(europa,SID==x)$Y)}))
                                   
#- Plot these rectangles + europa
plot(1,1,xlim=range(xrange),ylim=range(yrange),asp=1/lonLatRatio(mean(xrange),mean(yrange)),col="white",
     xlab="Longitude",ylab="Latitude")
plot(eurSP,col="darkgreen",add=T)
plot(rectsSP,add=T)
text(rects$SI_LONG,rects$SI_LATI,labels=rects$LE_RECT,cex=0.3)

#- Calculate proportion on land per ICES rectangle
propLand <- overlapPolygons(rectsSP,eurSP)
     
res <- expand.grid(1:480,1:1991)
res$total <- res[,1] * res[,2]