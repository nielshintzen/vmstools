
Grid2KML <- function(output.mat=output.mat, log.it=FALSE, what.quantity = 'effort',kmlfile='vms.kml',imagefile='vms.png')  {
#Takes the output from vmsGridCreate ie. when plotMap is set to FALSE.
#output.mat[["fishing"]] <- log(output.mat[["fishing"]])
dd <- output.mat@grid
d1 <- dd@cells.dim[1]
d2 <- dd@cells.dim[2]
fishing <- output.mat@data$fishing
mat <- matrix((fishing),byrow=TRUE,ncol=d1,nrow=d2)
mat <- t(mat[d2:1,])
if(log.it) {mat <- log(mat) }
bbox <- output.mat@bbox
xxx <- seq(bbox[1,1],bbox[1,2],length=d1)
yyy <- seq(bbox[2,1],bbox[2,2],length=d2)
rr <- range(mat[mat!='-Inf'],na.rm=T)
labs<-seq(rr[1],rr[2],length=9)
image(xxx,yyy,mat,zlim=c(rr[1],rr[2]),xlab="",ylab="",col=rainbow(9))
gd <- list(x=xxx,y=yyy,z=mat)
gd.1 <- as.SpatialGridDataFrame.im(as.im(gd))
proj4string(gd.1) <- CRS("+proj=longlat +datum=WGS84")
vms.kml <- GE_SpatialGrid(gd.1)
#tf <- tempfile(tmpdir=getwd())
png(file=imagefile, width=vms.kml$width, height=vms.kml$height, bg="transparent",res=576)
par(mar=c(0,0,0,0), xaxs="i", yaxs="i",cex=.25)
image(as.image.SpatialGridDataFrame(gd.1[1]), col=heat.colors(9),xlim=vms.kml$xlim, ylim=vms.kml$ylim)
kmlOverlay(vms.kml, kmlfile=kmlfile, imagefile=imagefile, name=what.quantity)

legend(x='bottomright', legend=as.character(round(labs)), pch = 22,border =heat.colors(length(labs)), pt.bg=heat.colors(length(labs)),
title=what.quantity, ncol=2,bg="white",pt.cex=.2 )

dev.off()
}


