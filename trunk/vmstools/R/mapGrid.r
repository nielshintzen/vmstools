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
                         , addLegend = TRUE
                         , legendx='bottomleft'
                         , legendncol = 1
                         , legendtitle = "fishing activity"
                         , plotPoints = FALSE
                         , outGridFile = ""  #name for output gridAscii
                         , outPlot = ""  #name for output png
                         , ... )
{

require(sp)
require(maptools)


xlim=c(we,ea)
ylim=c(so,no)

#dev.new()
breaks <- pretty(sGDF[[gridValName]],n=numCats)
#cols <- rev(heat.colors((length(breaks)-1)))
cols <- rev(terrain.colors((length(breaks)-1)))
image(sGDF, attr=gridValName, axes=FALSE,  col=cols, xlim=xlim, ylim=ylim,breaks=breaks)
library(mapdata)
map("worldHires", add=TRUE,col="darkgreen",fill=TRUE,bg="white")
box() # to put a box around the plot
#mtext(paste(gear,year),font=4,line=-1.5)

#legend(x='bottomleft', legend=breaks[1:(length(breaks)-1)], pch = 22, pt.bg=cols, title="fishing activity",bg="white",pt.cex=2 )
legend(x=legendx, legend=breaks[1:(length(breaks)-1)], pch = 22, pt.bg=cols, title=legendtitle, ncol=legendncol, bg="white",pt.cex=2 )

#to add points (can obscure grid)
if (plotPoints) points(sPDF, pch=".")

#to add plotTitle
if (plotTitle != "") mtext(plotTitle)



} #end of mapGrid