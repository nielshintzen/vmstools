nearest <- function(x,y){

  if(!class(x) %in% c("SpatialPoints","SpatialPointsDataFrame")) stop("No functionality build in for formats other than 'SpatialPoints' and 'SpatialPointsDataFrame'")
  if(!class(y) %in% c("SpatialLines","SpatialLinesDataFrame"))   stop("No functionality build in for formats other than 'SpatialLines' and 'SpatialLinesDataFrame'")

  #- Access lines
  lns         <- y@lines
  #- Access points
  pts         <- data.frame(coordinates(x))
  pts$ID      <- NA
  pts$idx     <- 1:nrow(pts); ptsOrig <- pts

  #- Turn lns into dataframe with coordinates only
  lncoords    <- data.frame(do.call(rbind,lapply(lns,function(z){do.call(rbind,coordinates(z))})))
  colnames(lncoords)[1:2] <- c("X1","X2")
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

  #- Setup ellipse calculation
  theta       <- seq(0, 2 * pi, len = 25)
  ID          <- numeric(nrow(pts)); ID[] <- NA
  IDln        <- lncoords$ID
  for(i in 1:nrow(pts)){
    ptsx      <- pts[i,1]; ptsy <- pts[i,2]
    llr       <- lonLatRatio(ptsx,ptsy)

    searchDegree <- 0.2; idx <- numeric()
    while(length(idx) == 0){
      searchDegree <- searchDegree + 0.2
      ell     <- cbind(cos(theta)*searchDegree +     ptsx,sin(theta)*searchDegree*llr + ptsy)
      ins     <- which(lncoords$X1 >= min(ell[,1]) & lncoords$X1 <= max(ell[,1]) &
                       lncoords$X2 >= min(ell[,2]) & lncoords$X2 <= max(ell[,2]))

      if(length(ins)>0){
        idx   <- point.in.polygon(lncoords$X1[ins],lncoords$X2[ins],ell[,1],ell[,2])
        idx   <- ins[which(idx == 1)]
      }
    }
    ID[i]     <- IDln[idx][which.min(distance(ptsx,ptsy,lncoords$X1[idx],lncoords$X2[idx]))]
  }
  ptsOrig$ID[pts$idx]      <- ID
return(ptsOrig)}


#- Example
#library(vmstools)
#load("D:/Repository/VMStools/polygons/BathymetryNorthSea/SpatialLinesDataFrame/bath5m.rdata")
#data(tacsat)
#
#y <- bath5m
#x <- SpatialPoints(data.frame(SI_LONG=tacsat$SI_LONG,SI_LATI=tacsat$SI_LATI))
#system.time(a <- nearest(x,y)) #37 minutes on  my machine
#
#tacsat$ID <- a$ID
#tacsat$DEPTH <- y@data[ac(tacsat$ID),"LLWS_M"]
#
#map("worldHires",xlim=c(0,8),ylim=c(50,56),fill=T, col="darkgreen");map.axes()
#plot(y,add=T,col="blue")
#for(st in 1:100){
#lines(coordinates(y@lines[[an(tacsat$ID[st])+1]])[[1]],col="red")
#print(paste(tacsat$ID[st],y@lines[[an(tacsat$ID[st])+1]]@ID))
#points(tacsat$SI_LONG[st],tacsat$SI_LATI[st],col="green",pch=19)
#Sys.sleep(0.5)
#}
#