pointInHarbour <- function(lon,lat,harbours,rowSize=30, returnNames=FALSE,saveHarbourList=TRUE){

    xharb     <- harbours$lon
    yharb     <- harbours$lat
    rharb     <- harbours$range
    harb      <- cbind(xharb,yharb,rharb)
    if("Description" %in% colnames(harbours)) rownames(harb) <- harbours$Description
    if("harbour" %in% colnames(harbours))     rownames(harb) <- harbours$harbour
    harb      <- orderBy(~xharb+yharb,data=harb)

    xys       <- data.frame(lon,lat)
    ordxys    <- order(xys$lon,xys$lat)
    lon       <- lon[ordxys]
    lat       <- lat[ordxys]

    nChunks   <- ceiling(length(lon)/rowSize)
     store     <- rep(0, length(lon))
    for(chunks in 1:nChunks){
      if(chunks == nChunks){
        x1    <- lon[(chunks*rowSize-rowSize+1):length(lon)]
        y1    <- lat[(chunks*rowSize-rowSize+1):length(lon)]
      } else {
          x1    <- lon[(chunks*rowSize-rowSize+1):(chunks*rowSize)]
          y1    <- lat[(chunks*rowSize-rowSize+1):(chunks*rowSize)]
        }

      xr        <- range(x1,na.rm=TRUE); xr <- c(xr[1]-0.05,xr[2]+0.05)
      yr        <- range(y1,na.rm=TRUE); yr <- c(yr[1]-0.05,yr[2]+0.05)
      res1      <- which(harb[,"xharb"] >= xr[1] & harb[,"xharb"] <= xr[2])
      res2      <- which(harb[,"yharb"] >= yr[1] & harb[,"yharb"] <= yr[2])
      res3      <- res1[which(is.na(pmatch(res1,res2))==FALSE)]

      if(length(res3)>0){
        for(hars in res3){
          #print(hars)
          x2  <- harb[hars,"xharb"]
          y2  <- harb[hars,"yharb"]

          pd  <- pi/180

          a1  <- sin(((y2-y1)*pd)/2)
          a2  <- cos(y1*pd)
          a3  <- cos(y2*pd)
          a4  <- sin(((x2-x1)*pd)/2)
          a   <- a1*a1+a2*a3*a4*a4

          c   <- 2*atan2(sqrt(a),sqrt(1-a));
          R   <- 6371;
          dx1 <- R*c

          res <- numeric(length(x1))
          idx <- which(dx1<=harb[hars,"rharb"])
          res[idx] <- 1
          if(returnNames){
            res[idx]  <- rownames(harb)[hars]  # overwrite '1' with the port names
            if(chunks==nChunks){
              idx2 <- idx[which(store[(chunks*rowSize-rowSize+1):length(lon)][idx] == "0")]
              store[(chunks*rowSize-rowSize+1):length(lon)][idx2] <- res[idx2]
            } else {
                idx2 <- idx[which(store[(chunks*rowSize-rowSize+1):(chunks*rowSize)][idx] == "0")]
                store[(chunks*rowSize-rowSize+1):(chunks*rowSize)][idx2] <- res[idx2]
              }
          } else {
              if(chunks==nChunks){
                store[(chunks*rowSize-rowSize+1):length(lon)]   <- store[(chunks*rowSize-rowSize+1):length(lon)]+res
              } else {
                  store[(chunks*rowSize-rowSize+1):(chunks*rowSize)] <- store[(chunks*rowSize-rowSize+1):(chunks*rowSize)]+res
                }
            }
        }
      }
    }
    if(returnNames == FALSE) store[which(store>0)] <- 1
    #Get order in tacsat back
    store[ordxys] <- store

    if(returnNames) store <- replace(store, store=="0", NA)
    if(saveHarbourList) write.table(harbours,file="harbourList_pointInHarbour.txt",append=FALSE,sep="\t")

return(store)}
    


