obs <- tacsat[round(seq(1,nrow(tacsat),length.out=30)),]
obs <- obs[,c(1,2,3,4,5,6)]
obs$OB_TYP <- "Benthos"
colnames(obs) <- c("OB_COU","OB_REF","SI_LATI","SI_LONG","SI_DATE","SI_TIME","OB_TYP")
obs$SI_TIME   <- ac(format(as.POSIXct(obs$SI_TIME,format="%H:%M")-runif(30,-60*20,60*20),format="%H:%M"))
obs$SI_LATI   <- jitter(obs$SI_LATI,factor=0.25)
obs$SI_LONG   <- jitter(obs$SI_LONG,factor=0.5)

a <- clipObs2Tacsat(tacsat,obs,method="grid",control.grid=list(resx=0.1,resy=0.05,gridBbox="obs"),temporalRange=c(-30,120),all.t=F)
a <- clipObs2Tacsat(tacsat,obs,method="euclidean",control.euclidean=list(threshold=1),temporalRange=c(-1e10,-1) ,all.t=F) 


clipObs2Tacsat <- function(tacsat,        #The tacsat dataset
                           obs,           #The observation dataset
                           method="grid", #The method used, on 'grid' or 'euclidean' distance
                           control.grid=list(spatGrid=NULL,resx=NULL,resy=NULL,gridBbox="obs"), #gridBbox: whether bounding box should come from tacsat or observations
                           control.euclidean=list(threshold=NULL), #all.t = all.tacsat
                           temporalRange=NULL,
                           all.t=F,
                           rowSize=1000
                           ){

tacsat      <- sortTacsat(tacsat)
obs         <- sortTacsat(obs)

  #- If you want to get the full tacsat dataset back, make a copy of it first
if(all.t){
  tacsat$ID   <- 1:nrow(tacsat)
  tacsatOrig  <- tacsat
}
if(!"SI_DATIM" %in% colnames(tacsat)) tacsat$SI_DATIM   <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME, sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
if(!"SI_DATIM" %in% colnames(obs))    obs$SI_DATIM      <- as.POSIXct(paste(obs$SI_DATE,    obs$SI_TIME,    sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")


#- Gridcell wanted, but not given yet, so create one
if(method == "grid" & is.null(control.grid$spatGrid) == T){
  if(is.null(control.grid$resx) == T | is.null(control.grid$resy) == T) stop("Method selected needs resx and resy statements")

  xrangeO      <- range(obs$SI_LONG,na.rm=T); xrangeT <- range(tacsat$SI_LONG,na.rm=T)
  yrangeO      <- range(obs$SI_LATI,na.rm=T); yrangeT <- range(tacsat$SI_LATI,na.rm=T)
  
  if(control.grid$gridBbox == "obs")     spatGrid    <- createGrid(xrangeO,yrangeO,control.grid$resx,control.grid$resx,type="SpatialGrid")
  if(control.grid$gridBbox == "tacsat")  spatGrid    <- createGrid(xrangeT,yrangeT,control.grid$resy,control.grid$resy,type="SpatialGrid")
}

#- Perform calculations on gridcell
if(method == "grid" & is.null(spatGrid) == F){
  sPDFObs     <- SpatialPointsDataFrame(data.frame(cbind(obs$SI_LONG,obs$SI_LATI)),data=obs)
  sPDFTac     <- SpatialPointsDataFrame(data.frame(cbind(tacsat$SI_LONG,tacsat$SI_LATI)),data=tacsat)
  resObs      <- overlay(spatGrid,sPDFObs)
  resTac      <- overlay(spatGrid,sPDFTac)

  idxObs      <- SpatialPoints(spatGrid)@coords[resObs,]
  idxTac      <- SpatialPoints(spatGrid)@coords[resTac,]
  
  obs$GR_LONG     <- idxObs[,1];  obs$GR_LATI     <- idxObs[,2]
  tacsat$GR_LONG  <- idxTac[,1];  tacsat$GR_LATI  <- idxTac[,2]

  tacsat$GR_ID<- resTac;
  obs$GR_ID   <- resObs;

  
  if(is.null(temporalRange)==F){

    res       <- do.call(c,lapply(as.list(1:nrow(obs)),function(x){     res        <- which(resTac == resObs[x]);
                                                                        restime    <- difftime(obs$SI_DATIM[x],tacsat$SI_DATIM[res],unit="mins");
                                                                        #retrn      <- ifelse(restime <= temporalRange[2] & restime >=temporalRange[1],resObs[x],NA)
                                                                        retrn      <- which(restime <= temporalRange[2] & restime >= temporalRange[1])
                                                       return(res[retrn])}))


    retrn       <- list(obs,tacsat[res,])
  } else {
      retrn       <- list(obs,tacsat)
    }
}

#- Perform calculation by Euclidian distance
if(method == "euclidean"){

  obs$GR_ID   <- 1:nrow(obs)

  #- Create storage of tacsat records to save
  totRes      <- cbind(numeric(), numeric())

  obsLon      <- obs$SI_LONG
  obsLat      <- obs$SI_LATI
  tacLon      <- tacsat$SI_LONG
  tacLat      <- tacsat$SI_LATI

  #- Chop it up into pieces to speed up the calculations
  nChunkObs   <- ceiling(nrow(obs)    /rowSize)
  nChunkTac   <- ceiling(nrow(tacsat) /rowSize)
  for(iNO in 1:nChunkObs){
    if(iNO == nChunkObs){
      ox <- obsLon[(iNO*rowSize - rowSize + 1):length(obsLon)]
      oy <- obsLat[(iNO*rowSize - rowSize + 1):length(obsLat)]
    } else {
        ox <- obsLon[(iNO*rowSize - rowSize + 1):(iNO*rowSize)]
        oy <- obsLat[(iNO*rowSize - rowSize + 1):(iNO*rowSize)]
      }
    for(iNT in 1:nChunkTac){
      if(iNT == nChunkTac){
        tx <- tacLon[(iNT*rowSize - rowSize + 1):length(tacLon)]
        ty <- tacLat[(iNT*rowSize - rowSize + 1):length(tacLat)]
      } else {
          tx <- tacLon[(iNT*rowSize - rowSize + 1):(iNT*rowSize)]
          ty <- tacLat[(iNT*rowSize - rowSize + 1):(iNT*rowSize)]
        }

      #- Check if the length of both sets are equal or not
      if(iNO == nChunkObs | iNT == nChunkTac){

        #- Get both the minimum distance, but also the index of the tacsat record associated
        res                                                   <- do.call(rbind,
                                                                         lapply(as.list(1:length((iNO*rowSize - rowSize +1):ifelse(iNO == nChunkObs,length(obsLon),(iNO*rowSize)))),
                                                                                function(x){
                                                                                    #- Get the row numbers of the full observation and tacsat set used here
                                                                                    obsRows     <- (iNO*rowSize - rowSize +1):ifelse(iNO == nChunkObs,length(obsLon),(iNO*rowSize))
                                                                                    if(length(tx)<rowSize){  tacRows <- ((iNT*rowSize - rowSize + 1):(length(tacLon))) } else { tacRows <- ((iNT*rowSize - rowSize + 1):(iNT*rowSize))}

                                                                                    #- Calculate the distance between the observation and tacsat dataset
                                                                                    distObsTac  <- distance(ox[x],oy[x],tx,ty);
                                                                                    if(is.null(control.euclidean$threshold)==F){ idx <- which(distObsTac <= control.euclidean$threshold)} else { idx <- 1:length(distObsTac) }

                                                                                    restime     <- difftime(obs$SI_DATIM[obsRows[x]],tacsat$SI_DATIM[tacRows[idx]],unit="mins")
                                                                                    if(is.null(temporalRange)==F){ retrn       <- which(restime <= temporalRange[2] & restime >= temporalRange[1])
                                                                                    } else { retrn <- 1:length(restime) }
                                                                                    if(length(tacRows[idx[retrn]])>0){ toReturn <- cbind(tacRows[idx[retrn]],obs$GR_ID[obsRows[x]])
                                                                                    } else { toReturn <- cbind(NA,NA) }
                                                                                    
                                                                                return(toReturn)}))



        totRes  <- rbind(totRes,na.omit(res))
      }
      if(iNO != nChunkObs & iNT != nChunkTac){

       obsRows            <- ((iNO*rowSize - rowSize +1):(iNO*rowSize))
       tacRows            <- (((iNT*rowSize - rowSize + 1):(iNT*rowSize)))
       res                <- outer(1:rowSize,1:rowSize,FUN=
                                   function(x,y){
                                      distObsTac  <- distance(ox[x],oy[x],tx[y],ty[y])
                                   return(distObsTac)})
                                
       idx              <- apply(res,1,function(x){return(x <= control.euclidean$threshold)})
       idx              <- which(idx == T,arr.ind=T)
       restime          <- difftime(obs$SI_DATIM[obsRows[idx[,1]]],tacsat$SI_DATIM[tacRows[idx[,2]]],unit="mins")

       if(is.null(temporalRange)==F){ retrn       <- which(restime <= temporalRange[2] & restime >= temporalRange[1])
       } else { retrn <- 1:length(restime) }
       res              <- cbind(tacRows[idx[retrn,2]],obs$GR_ID[obsRows[idx[,1]]])
       totRes           <- rbind(totRes,na.omit(res))
      }
      
      
    }#End iNT loop
  }#End iNO loop
  tacsat$GR_ID <- NA
  tacsat$GR_ID[totRes[,1]] <-totRes[,2]


  retrn   <- list(obs,tacsat[totRes[,1],])
}#End method euclidean

if(all.t) retrn[[2]] <- merge(retrn[[2]],tacsatOrig,by=colnames(tacsatOrig),all=T)

return(retrn)}
