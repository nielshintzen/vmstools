#' Find points in harbour within specified range
#' 
#' Method to find the gps positions given with tacsat data that are situated
#' within a range of a port.
#' 
#' The method returns the index of points that are within a harbour area, given
#' the midpoints of the harbours and a range (in km) from these midpoints.
#' 
#' @param lon Longitudinal positions of the TACSAT formatted data
#' @param lat Latitudinal positions of teh TACSAT formatted data
#' @param harbour Latitudinal and Longitudinal position of the harbour and
#' outer range from midpoint of harbour
#' @param returnNames Logical: return the name of the harbour instead of 1 / 0
#' indicating if it is inside or outside the harbour. Default to FALSE
#' @param saveHarbourList Logical: writing harbour list used in function to
#' file. Default to TRUE
#' @author Niels T. Hintzen
#' @seealso \code{\link{distance}}, \code{\link{lonLatRatio}},
#' \code{\link{sortTacsat}}, \code{\link{filterTacsat}},
#' \code{\link{mergeEflalo2Tacsat}}
#' @references EU lot 2 project
#' @examples
#' 
#' data(eflalo)
#' data(tacsat)
#' data(euharbours); euharbours <- harbours
#' 
#' #-Remove duplicated records from tacsat
#' myf       <- paste(tacsat$VE_REF,tacsat$SI_LATI,tacsat$SI_LONG,
#'                    tacsat$SI_DATE,tacsat$SI_TIME);
#' tacsat    <- tacsat[!duplicated(myf),];
#' 
#' 
#' #-Find all the gps locations that are located within the port area
#' idx <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,
#'                       harbours=harbours,returnNames=TRUE)
#' print(head(idx))
#' getwd() #in this directory, the harbour list will be written to disk
#' idx <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,
#'                       harbours=harbours,saveHarbourList=TRUE)
#' 
#' idx <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,
#'                       harbours=harbours)
#' idx <- which(idx==1)
#' 
#' #-Plot these port locations on a map
#' library(maps); library(mapdata)
#' #map the world, but plot only the northsea by lon and lat limits,
#' # in high resolution
#' xrange <- range(tacsat$SI_LONG[idx])
#' yrange <- range(tacsat$SI_LATI[idx])
#' 
#' map('worldHires',xlim=xrange,ylim=yrange,col="darkgreen",fill=TRUE,
#'     resolution=1, bg="white", border=0)
#' map.axes();
#' mtext("Longitude",1,line=-2,outer=TRUE,cex=1.2,font=2)
#' mtext("Latitude",2,line=-2,outer=TRUE,cex=1.2,font=2)
#' 
#' points(tacsat$SI_LONG[idx],tacsat$SI_LATI[idx],cex=0.1,pch=16,col="red")
#' 
#' @export pointInHarbour
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
    


