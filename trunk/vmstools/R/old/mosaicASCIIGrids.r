## mosaicASCIIGrids.r
## by Fabrizio Manco, 10/11/2010
## merge a list of ASCII grids

mosaicASCIIGrids <- function (ascgridfilelist,
                             outputFileName)

{
  require(raster)
  require(maptools)
  grid1<-readAsciiGrid(ascgridfilelist[1])
  mergedRaster<-raster(grid1)
  rm(grid1)
  
  for (x in 2:length(ascgridfilelist)){
    grid2<-readAsciiGrid(ascgridfilelist[x])
    raster2<-raster(grid2)
    
    # do the merging and export as a ascii file
    mergedRaster<-mosaic(mergedRaster, raster2, fun=sum,tolerance=0.1)
    }
  writeRaster(mergedRaster, outputFileName, format="ascii", overwrite=TRUE)
}


