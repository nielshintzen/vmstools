#mapGrid.r
#andy south 12/2/09

#to map grids input as SGDF
#** OR  maybe just provide a wrapper to mapGriddedData instead ?
#would need to add a better worldmap into that

mapGrid <- function( sGDF
                         , sPDF
                         , we=""
                         , ea=""
                         , so=""
                         , no=""
                         , gridValName="fishing"
                         , plotTitle = ""
                         , numCats = 5
                         , paletteCats = "heat.colors"
                         , addLegend = TRUE
                         , legendx='bottomleft'
                         , legendncol = 1
                         , legendtitle = "fishing activity"
                         , plotPoints = FALSE
                         , colPoints =1
                         , colLand = 'sienna'
                         , addScale = TRUE
                         , outGridFile = ""  #name for output gridAscii
                         , outPlot = ""  #name for output png
                         , ... )
{

require(sp)
require(maptools)

par(mar=c(4,6,1,1))
 
xlim0=c(we,ea)
ylim0=c(so,no)

lstargs <- list(...)

#dev.new()
if(length(lstargs$breaks0)==0) {
      breaks0 <- pretty(sGDF[[gridValName]],n=numCats)
      } else{
      breaks0 <- lstargs$breaks0
      }

# rainbow, heat.colors, etc.
cols <- rev(do.call(paletteCats, list(length(breaks0)-1)))


library(mapdata)
map("worldHires", add=FALSE,col=colLand,fill=TRUE, bg="white",  xlim=xlim0 + c(+0.1,-0.1), ylim=ylim0 + c(+0.1,-0.1), 
regions=c('uk','ireland','france','germany','netherlands', 'norway','belgium',
'spain','luxembourg','denmark', 'sweden','iceland', 'portugal','italy','sicily','ussr','sardinia','albania','monaco','turkey','austria',
'switzerland','czechoslovakia','finland','libya', 'hungary','yugoslavia','poland','greece','romania','bulgaria', 'slovakia','morocco',
'tunisia','algeria','egypt' ))

  im <-  as.image.SpatialGridDataFrame(sGDF, attr=2)
   image(im$x,im$y,im$z,  axes=FALSE, col=cols,  breaks = breaks0, 
   xlim=xlim0 , ylim=ylim0, add=TRUE   )
   

#map("worldHires", add=TRUE, col=colLand, fill=TRUE, bg="white",  xlim=xlim0 + c(+1,-1), ylim=ylim0 + c(+1,-1), 
map("worldHires", add=TRUE, col=colLand, fill=TRUE, bg="white",  xlim=xlim0 , ylim=ylim0 , 
regions=c('uk','ireland','france','germany','netherlands', 'norway','belgium',
'spain','luxembourg','denmark', 'sweden','iceland', 'portugal','italy','sicily','ussr','sardinia','albania','monaco','turkey','austria',
'switzerland','czechoslovakia','finland','libya', 'hungary','yugoslavia','poland','greece','romania','bulgaria', 'slovakia','morocco',
'tunisia','algeria','egypt' ),   mar=c(6,6,2,2))


box() # to put a box around the plot
#mtext(paste(gear,year),font=4,line=-1.5)
axis(1)
axis(2, las=2)
 if(we>0){
   mtext("Degree East", side=1, line=2)
   } else{
   mtext("Degree West", side=1, line=2)
   }
 if(no>0){
   mtext("Degree North", side=2, line=3)
 } else{
   mtext("Degree South", side=2, line=3)
   }

 

# add a scale
if(addScale) map.scale(x=xlim0[2]-(xlim0[2]-xlim0[1])/2, y=ylim0[1], ratio=FALSE)

#to add points (can obscure grid)
if (plotPoints) points(sPDF, pch=".", col=colPoints)

#legend(x='bottomleft', legend=breaks[1:(length(breaks)-1)], pch = 22, pt.bg=cols, title="fishing activity",bg="white",pt.cex=2 )
legend(x=legendx, legend=breaks0[1:(length(breaks0)-1)], pch = 22, pt.bg=cols, title=legendtitle, ncol=legendncol, bg="white",pt.cex=2 )


#to add plotTitle
if (plotTitle != "") mtext(plotTitle)



} #end of mapGrid