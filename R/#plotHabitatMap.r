plotHabitatMap <- function(read.in.shapefile=TRUE,what.habitat.grouping="Grouped",min.lon=-5, max.lon=10, min.lat=48, max.lat=62,OS='Windoze')
{

require(maps)
require(mapdata)
#what.habitat.grouping = "substrate"
#what.habitat.grouping = "SubsGroups"
#what.habitat.grouping ="EUSMRegion"
#what.habitat.grouping ="SubsGrpPlu"

if(read.in.shapefile){
if(OS == "Windoze") habitat   <- readShapePoly("n:/Projecten/VMS tools/10 Papers/VMS impact indicators/data/HabitatMaps/ModelledSeabedHabitats.shp")
if(OS == "Linux") habitat   <- readShapePoly("/media/n/Projecten/VMS tools/10 Papers/VMS impact indicators/data/HabitatMaps/ModelledSeabedHabitats.shp")
}


print(names(habitat)) # print out potential polygons to examine
print(unique(habitat[[what.habitat.grouping]]))
#colors()[grep('brown',colors())]

map("worldHires",xlim=c(min.lon,max.lon),ylim=c(min.lat,max.lat),fill=TRUE,col="darkgreen"); map.axes()

sp <- SpatialPolygons(habitat@polygons)

#sb <- grep("Seabed", habitat[[what.habitat.grouping]])
#msm <- grep("Mud to sandy mud", habitat[[what.habitat.grouping]])
#sms <- grep("Sand to muddy sand", habitat[[what.habitat.grouping]])
#rock <- grep("Rock or other hard substrata", habitat[[what.habitat.grouping]])
#
#plot(sp[sms],add=TRUE,col="yellow",border="yellow")
#plot(sp[rock],add=TRUE,col="green",border="green")
#
#plot(sp[sb],add=TRUE,col="thistle",border="thistle")
#plot(sp[msm],add=TRUE,col="brown",border="brown")
#

if (what.habitat.grouping == "Grouped")  {
colrs <- c(
"lightblue","tan","yellow","lightgray","green","red",
"blue","saddlebrown","yellow1","lightslategray","green1",
"blue1","sandybrown","yellow2","slategrey","green2",
"darkblue","brown4","lightyellow4","gray1","darkolivegreen")
 }

if (what.habitat.grouping == "substrate")  {
colrs <- c(
"lightblue","brown","yellow","slategray","green","red",
"blue")
 }

if (what.habitat.grouping == "SubsGrpPlu")  {
colrs <- c(
"lightblue","brown","yellow","green")
 }


lcolrs <- length(colrs)+1 
print(colrs)


plot(sp,col=colrs,add=TRUE,border=colrs) 

map("worldHires",xlim=c(min.lon,max.lon),ylim=c(min.lat,max.lat),fill=TRUE,add=TRUE,col="darkgreen")
legend("bottomleft",legend=c("Land",ac(unique(habitat[[what.habitat.grouping]]))),

        lwd=c(1,1,1,1),lty=rep(NA,10),pch=rep(22,lcolrs),ncol=2,
        col=rep("black",lcolrs),
        pt.bg=c("darkgreen",colrs),pt.cex=rep(3,lcolrs),box.lty=1,box.lwd=2,
        y.intersp=1.5,title="Substrate type",bg="white",cex=0.5)
        
        }


#OS <- 'Linux'
#plotHabitatMap(read.in.shapefile=TRUE,what.habitat.grouping="SubsGrpPlu",min.lon= -5,max.lon=10, min.lat=45,max.lat=65)


map("worldHires",xlim=c(-5,10),ylim=c(48,62),fill=TRUE,col="darkgreen");
map.axes()

sp <- SpatialPolygons(habitat@polygons)
plot(sp,col=colrs[an(habitat$substrate)],
     border=colrs[an(habitat$substrate)],add=TRUE)
     
legend("bottomleft",legend=c("Land",levels(habitat$substrate)),
        lwd=c(1,1,1,1),lty=rep(NA,10),pch=rep(22,lcolrs),ncol=2,
        col=rep("black",lcolrs), pt.bg=c("darkgreen",colrs),pt.cex=rep(3,lcolrs),box.lty=1,box.lwd=2,
        y.intersp=1.5,title="Substrate type",bg="white",cex=0.5)

for(i in levels(habitat$substrate)){
x11()
        
idx <- which(habitat$substrate == i)
sp2 <- SpatialPolygons(habitat@polygons[idx])
map("worldHires",xlim=c(-5,10),ylim=c(48,62),fill=TRUE,col="darkgreen");
title(i)
map.axes()
plot(sp2,col=colrs[an(habitat$substrate[idx])[order(sp2@plotOrder)]],
     border=colrs[an(habitat$substrate[idx])[order(sp2@plotOrder)]],add=TRUE)
}


pnt <- c(4.8,54.2) #Mud to sandy mut
pnt <- c(0.7,58.72) #Mixed sediment
pnt <- c(4.5,59.3) #Till
pnt <- c(2,56) #Sand to muddy sand
insidePols <- rep(0,length(habitat@polygons))
for(i in 1:length(habitat@polygons)){
  for(j in 1:length(habitat@polygons[[i]]@Polygons)){
    pol <- coordinates(habitat@polygons[[i]]@Polygons[[j]])
    insidePols[i] <- max(insidePols[i],point.in.polygon(pnt[1],pnt[2],pol[,1],pol[,2]))
  }
}
idx <- which(insidePols == 1)
habitat$substrate[idx]
habitat$substrate[idx][which.max(order(habitat@plotOrder)[idx])]