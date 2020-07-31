nearest <- function(x,y,rowSize=30){
  require(fields)

  if(!class(x) %in% c("SpatialPoints","SpatialPointsDataFrame")) stop("No functionality build in for formats other than 'SpatialPoints' and 'SpatialPointsDataFrame'")
  if(!class(y) %in% c("SpatialLines","SpatialLinesDataFrame"))   stop("No functionality build in for formats other than 'SpatialLines' and 'SpatialLinesDataFrame'")

  if(length(unique(lapply(y@lines,function(z){z@ID}))) != length(y)) stop("Each line segment does not contain a unique ID")

  #- Access lines
  lns         <- y@lines
  #- Access points
  pts         <- data.frame(coordinates(x))
  pts$ID      <- NA
  pts$idx     <- 1:nrow(pts); ptsOrig <- pts

  #- Turn lns into dataframe with coordinates only
  lncoords    <- data.frame(do.call(rbind,lapply(lns,function(z){do.call(rbind,coordinates(z))})))
  IDs         <- unlist(lapply(lns,function(z){z@ID}))
  lnlengths   <- unlist(lapply(lapply(lns,function(z){do.call(rbind,coordinates(z))}),nrow))
  lncoords$ID <- rep(IDs,lnlengths)

  #- Remove already points with no overlap with lines object
  xbox        <- bbox(x)
  ybox        <- bbox(y)

  remx <- which(pts[,1] < ybox[1,1] | pts[,1] > ybox[1,2] | pts[,2] < ybox[2,1] | pts[,2] > ybox[2,2])

  if(length(remx)>0) pts      <- pts[-remx,]

  #- Calculate nearest distance of each point to line
  lncoords    <- orderBy(~X1+X2,data=lncoords)
  pts         <- orderBy(~SI_LONG+SI_LATI,data=pts)

  out         <  - rdist.earth ( ozone2$lon.lat)


  #- Setup ellipse calculation
  theta       <- seq(0, 2 * pi, len = 25)
  for(i in 1:nrow(pts)){

  start.time <- Sys.time()
  for(i in 1:10000){
    llr       <- lonLatRatio(pts[i,1],pts[i,2])

    searchDegree <- 0.025; idx <- numeric()
    while(length(idx) == 0){
      searchDegree <- searchDegree + 0.025
      ell       <- cbind(cos(theta)*searchDegree +     pts[i,1],
                         sin(theta)*searchDegree*llr + pts[i,2])
      xr        <- c(min(ell[,1]),max(ell[,1])); yr <- c(min(ell[,2]),max(ell[,2]))
      ins       <- which(lncoords[,1] >= xr[1] & lncoords[,1] <= xr[2] &
                         lncoords[,2] >= yr[1] & lncoords[,2] <= yr[2])

      if(length(ins)>0){
        idx       <- point.in.polygon(lncoords[ins,1],lncoords[ins,2],ell[,1],ell[,2])
        idx       <- ins[which(idx == 1)]
      }
    }
    pts$ID[i]   <- lncoords$ID[idx][which.min(distance(pts[i,1],pts[i,2],lncoords[idx,1],lncoords[idx,2]))]
  }
  difftime(Sys.time(),start.time)



return(pts)}


#- Example
library(vmstools)
load("D:/Repository/VMStools/polygons/BathymetryNorthSea/SpatialLinesDataFrame/bath5m.rdata")
data(tacsat)

y <- bath5m
x <- SpatialPoints(data.frame(SI_LONG=tacsat$SI_LONG,SI_LATI=tacsat$SI_LATI))

st <- 8610
map("worldHires",xlim=c(0,8),ylim=c(50,56),fill=T, col="darkgreen");map.axes()
plot(y,add=T,col="blue")
lines(coordinates(y@lines[[an(pts$ID[st])+1]])[[1]],col="red")
print(paste(pts$ID[st],y@lines[[an(pts$ID[st])+1]]@ID))
points(pts[st,1],pts[st,2],col="green",pch=19)

      ell <- cbind(cos(theta)*(xr[2]-sum(xr)/2)*1/cos(0.25*pi) + sum(xr)/2,
                   sin(theta)*(yr[2]-sum(yr)/2)*1/sin(0.75*pi) + sum(yr)/2)


      plot(ell,type="l")
      points(xr,yr,pch=19,col="red")
      points(subPts)
      
      
nearest <- function(x,y,rowSize=30){

  if(!class(x) %in% c("SpatialPoints","SpatialPointsDataFrame")) stop("No functionality build in for formats other than 'SpatialPoints' and 'SpatialPointsDataFrame'")
  if(!class(y) %in% c("SpatialLines","SpatialLinesDataFrame"))   stop("No functionality build in for formats other than 'SpatialLines' and 'SpatialLinesDataFrame'")

  #- Access lines
  lns         <- y@lines
  #- Access points
  pts         <- data.frame(coordinates(x))
  pts$ID      <- NA
  pts$idx     <- 1:nrow(pts); ptsOrig <- pts
  pts$minDist <- 1e10

  #- Turn lns into dataframe with coordinates only
  lncoords    <- data.frame(do.call(rbind,lapply(lns,function(z){do.call(rbind,coordinates(z))})))
  IDs         <- unlist(lapply(lns,function(z){z@ID}))
  lnlengths   <- unlist(lapply(lapply(lns,function(z){do.call(rbind,coordinates(z))}),nrow))
  lncoords$ID <- rep(IDs,lnlengths)

  #- Remove already points with no overlap with lines object
  xbox        <- bbox(x)
  ybox        <- bbox(y)

  remx <- which(pts[,1] < ybox[1,1] | pts[,1] > ybox[1,2] | pts[,2] < ybox[2,1] | pts[,2] > ybox[2,2])
  remy <- which(lncoords[,1] < xbox[1,1] | lncoords[,1] > xbox[1,2] | lncoords[,2] < xbox[2,1] | lncoords[,2] > xbox[2,2])

  if(length(remy)>0) lncoords <- lncoords[-remy,]
  if(length(remx)>0) pts      <- pts[-remx,]

  #- Calculate nearest distance of each point to line
  lncoords    <- orderBy(~X1+X2,data=lncoords)
  pts         <- orderBy(~SI_LONG+SI_LATI,data=pts)




  coslat1 <- cos((pts[, 2] * pi)/180)
  sinlat1 <- sin((pts[, 2] * pi)/180)
  coslon1 <- cos((pts[, 1] * pi)/180)
  sinlon1 <- sin((pts[, 1] * pi)/180)

  coslat2 <- cos((lncoords[, 2] * pi)/180)
  sinlat2 <- sin((lncoords[, 2] * pi)/180)
  coslon2 <- cos((lncoords[, 1] * pi)/180)
  sinlon2 <- sin((lncoords[, 1] * pi)/180)


  nChunkPts <- ceiling(nrow(pts)/rowSize)
  nChunkLin <- ceiling(nrow(lncoords)/rowSize)
  minDist   <- numeric(nrow(pts)); minDist[] <- 1e10
  ID        <- numeric(nrow(pts)); ID <- NA
  IDln      <- lncoords$ID

start.time <- Sys.time()
for(i in 1:nChunkPts){
  print(i)
  for(j in 1:nChunkLin){

    if(i == nChunkPts & j != nChunkPts){
    rowPts <- (i * rowSize - rowSize + 1):nrow(pts)
    rowLns <- (j * rowSize - rowSize + 1):(j * rowSize)
    }
    if(i == nChunkPts & j == nChunkPts){
     rowPts <- (i * rowSize - rowSize + 1):nrow(pts)
     rowLns <- (j * rowSize - rowSize + 1):nrow(lncoords)
    }
    if(i != nChunkPts & j == nChunkPts){
     rowPts <- (i * rowSize - rowSize + 1):(i * rowSize)
     rowLns <- (j * rowSize - rowSize + 1):nrow(lncoords)
    }
    if(i != nChunkPts & j != nChunkPts){
     rowPts <- (i * rowSize - rowSize + 1):(i * rowSize)
     rowLns <- (j * rowSize - rowSize + 1):((j * rowSize))
    }

    pp <-   cbind(coslat1[rowPts] * coslon1[rowPts], coslat1[rowPts] * sinlon1[rowPts], sinlat1[rowPts]) %*%
          t(cbind(coslat2[rowLns] * coslon2[rowLns], coslat2[rowLns] * sinlon2[rowLns], sinlat2[rowLns]))
    pp[abs(pp)>1]             <- -1
    dists                     <- 6371 * acos(pp)
    tmpminDist                <- apply(dists,1,min)
    wminDist                  <- apply(dists,1,which.min)
    idx                       <- which(tmpminDist < minDist[rowPts])
    minDist[rowPts[idx]]      <- tmpminDist[idx]
    ID[rowPts[idx]]           <- IDln[rowLns][wminDist[idx]]
   }
}
difftime(Sys.time(),start.time)





  #- Setup ellipse calculation
  theta       <- seq(0, 2 * pi, len = 25)
  for(i in 1:nrow(pts)){

  start.time <- Sys.time()
  for(i in 1:10000){
    llr       <- lonLatRatio(pts[i,1],pts[i,2])

    searchDegree <- 1; idx <- numeric()
    while(length(idx) == 0){
      searchDegree <- searchDegree + 0.025
      ell       <- cbind(cos(theta)*searchDegree +     pts[i,1],
                         sin(theta)*searchDegree*llr + pts[i,2])
#      ins       <- which(lncoords[,1] >= (pts[i,1]-2*(searchDegree-1)) & lncoords[,1] <= (pts[i,1]+2*(searchDegree-1)) &
#                         lncoords[,2] >= (pts[i,2]-2*(searchDegree-1)) & lncoords[,2] <= (pts[i,2]+2*(searchDegree-1)))
#      if(length(ins)>0){
#        idx       <- point.in.polygon(lncoords[ins,1],lncoords[ins,2],ell[,1],ell[,2])
        idx       <- point.in.polygon(lncoords[,1],lncoords[,2],ell[,1],ell[,2])
        idx       <- ins[which(idx == 1)]
      }
    }
    pts$ID[i]   <- lncoords$ID[idx][which.min(distance(pts[i,1],pts[i,2],lncoords[idx,1],lncoords[idx,2]))]
  }
  difftime(Sys.time(),start.time)


  #- Look over checking which line is closest in blocks
  nChunks <- ceiling(nrow(pts)/rowSize)
  for (chunks in 1:nChunks){
    if (chunks == nChunks) {
      rows <- (chunks * rowSize - rowSize + 1):nrow(pts)
    } else {
        rows <- (chunks * rowSize - rowSize + 1):(chunks * rowSize)
      }
    subPts <- pts[rows,]

    a<- rdist.earth(subPts[,1:2],lncoords[,1:2])



    idx <- numeric(); searchDegree <- 0.025
    while(length(idx) == 0){
      searchDegree <- searchDegree + 0.025
      xr  <- range(subPts[,1], na.rm = TRUE)
      xr  <- c(xr[1] - searchDegree, xr[2] + searchDegree)
      yr  <- range(subPts[,2], na.rm = TRUE)
      llr <- lonLatRatio(sum(xr)/2,sum(yr)/2) #Lon-Lat ratio
      yr  <- c(yr[1] - searchDegree*llr, yr[2] + searchDegree*llr)

      idx <- which(lncoords[,1] >= xr[1] & lncoords[,1] <= xr[2] &
                   lncoords[,2] >= yr[1] & lncoords[,2] <= yr[2])
    }

    dists <- outer(1:length(rows),1:length(idx),function(i,j){
              distance(subPts[i,1],subPts[i,2],lncoords[idx,][j,1],lncoords[idx,][j,2])}))
    if(chunks ==  nChunks){
      pts$ID[rows] <- lncoords$ID[idx][apply(dists,1,which.min)]
    } else {
        pts$ID[rows] <- lncoords$ID[idx][apply(dists,1,which.min)]
      }
  }
  ptsOrig$ID[pts$idx] <- pts$ID

return(pts)}


#- Example
library(vmstools)
load("D:/Repository/VMStools/polygons/BathymetryNorthSea/SpatialLinesDataFrame/bath5m.rdata")
data(tacsat)

y <- bath5m
x <- SpatialPoints(data.frame(SI_LONG=tacsat$SI_LONG,SI_LATI=tacsat$SI_LATI))

st <- 8610
map("worldHires",xlim=c(0,8),ylim=c(50,56),fill=T, col="darkgreen");map.axes()
plot(y,add=T,col="blue")
lines(coordinates(y@lines[[an(pts$ID[st])+1]])[[1]],col="red")
print(paste(pts$ID[st],y@lines[[an(pts$ID[st])+1]]@ID))
points(pts[st,1],pts[st,2],col="green",pch=19)

      ell <- cbind(cos(theta)*(xr[2]-sum(xr)/2)*1/cos(0.25*pi) + sum(xr)/2,
                   sin(theta)*(yr[2]-sum(yr)/2)*1/sin(0.75*pi) + sum(yr)/2)


      plot(ell,type="l")
      points(xr,yr,pch=19,col="red")
      points(subPts)