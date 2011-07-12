plot.habitat.map <- function(read.in.shapefile=TRUE,what.habitat.grouping="Grouped",min.lon=-5, max.lon=10, min.lat=48, max.lat=62,OS='Windoze')
{
#what.habitat.grouping = "substrate"
#what.habitat.grouping = "SubsGroups"
#what.habitat.grouping ="EUSMRegion"
#what.habitat.grouping ="SubsGrpPlu"

if(read.in.shapefile){
if(OS == "Windoze") habitat   <- readShapePoly("n:/Projecten/VMS tools/10 Papers/VMS impact indicators/data/ModelledSeabedHabitats.shp")
if(OS == "Linux") habitat   <- readShapePoly("/media/n/Projecten/VMS tools/10 Papers/VMS impact indicators/data/ModelledSeabedHabitats.shp")
}




print(names(habitat)) # print out potential polygons to examine
print(unique(habitat[[what.habitat.grouping]]))
#colors()[grep('brown',colors())]

map("worldHires",xlim=c(min.lon,max.lon),ylim=c(min.lat,max.lat),fill=T,col="darkgreen"); map.axes()

sp <- SpatialPolygons(habitat@polygons)

#sb <- grep("Seabed", habitat[[what.habitat.grouping]])
#msm <- grep("Mud to sandy mud", habitat[[what.habitat.grouping]])
#sms <- grep("Sand to muddy sand", habitat[[what.habitat.grouping]])
#rock <- grep("Rock or other hard substrata", habitat[[what.habitat.grouping]])
#
#plot(sp[sms],add=T,col="yellow",border="yellow")
#plot(sp[rock],add=T,col="green",border="green")
#
#plot(sp[sb],add=T,col="thistle",border="thistle")
#plot(sp[msm],add=T,col="brown",border="brown")
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


plot(sp,col=colrs,add=T,border=colrs) 

map("worldHires",xlim=c(min.lon,max.lon),ylim=c(min.lat,max.lat),fill=T,add=T,col="darkgreen")
legend("bottomleft",legend=c("Land",ac(unique(habitat[[what.habitat.grouping]]))),

        lwd=c(1,1,1,1),lty=rep(NA,10),pch=rep(22,lcolrs),ncol=2,
        col=rep("black",lcolrs),
        pt.bg=c("darkgreen",colrs),pt.cex=rep(3,lcolrs),box.lty=1,box.lwd=2,
        y.intersp=1.5,title="Substrate type",bg="white",cex=0.5)
        
        }



plot.habitat.map(read.in.shapefile=FALSE,what.habitat.grouping="SubsGrpPlu",min.lon= -5,max.lon=10, min.lat=45,max.lat=65)
