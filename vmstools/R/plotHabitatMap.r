plot.habitat.map <- function(min.lon=-5, max.lon=10, min.lat=48, max.lat=62,OS='Windoze')
{

if(OS == "Windoze") habitat   <- readShapePoly("n:/Projecten/VMS tools/10 Papers/VMS impact indicators/data/ModelledSeabedHabitats.shp")
if(OS == "Linux") habitat   <- readShapePoly("/media/n/Projecten/VMS tools/10 Papers/VMS impact indicators/data/ModelledSeabedHabitats.shp")

map("worldHires",xlim=c(min.lon,max.lon),ylim=c(min.lat,max.lat),fill=T,col="darkgreen"); map.axes()
plot(SpatialPolygons(habitat@polygons),col=an(habitat$substrate)+1,add=T,border=0)
map("worldHires",xlim=c(-5,10),ylim=c(48,62),fill=T,col="darkgreen",add=T); map.axes()
legend("bottomright",legend=c("Land",ac(unique(habitat$substrate))),
        lwd=c(1,1,1,1),lty=rep(NA,10),pch=rep(22,8),ncol=1,
        col=rep("black",8),
        pt.bg=c("darkgreen",an(unique(habitat$substrate))+1),pt.cex=rep(3,8),box.lty=1,box.lwd=2,
        y.intersp=1.2,title="Substrate type",bg="white",cex=0.8)

}

# plot.habitat.map()
