nestedGrid <- function(tacsat,resx,resy,maxstep = 10, n = 20,control=list(clm=NULL,FUN=NULL)){
  lon       <- tacsat$SI_LONG; lat <- tacsat$SI_LATI
  poly_wkt  <- rep(NA,length(lon))
  idx       <- 1:length(lon)
  i         <- 1
  while(i <= maxstep & length(idx) > 0){
    x       <- 2^floor((i-1)/2)
    y       <- 2^floor((i)/2)
    poly    <- polyDef(lon[idx], lat[idx], resx/x, resy/y)
    count   <- table(poly)

    idxin   <- which(poly  %in% names(count[which(count >= n)]))
    idxout  <- which(poly  %in% names(count[which(count < n)]))
    poly_wkt[idx[idxout]] <- poly[idxout]

    #- Final rounds
    if(i == maxstep){
      idxout  <- 1:length(poly)
      poly_wkt[idx[idxout]] <- poly[idxout]
    }

    idx       <- idx[idxin]
    i <- i + 1
  }

  #- Create SpatialPolygons
  uniquePols  <- poly_wkt[!duplicated(poly_wkt)]
  Pol         <- lapply(as.list(1:length(uniquePols)),
                  function(i){eval(parse(text=uniquePols[i]))})
  Pols        <- lapply(as.list(1:length(uniquePols)),function(x){Polygons(list(Pol[[x]]),ID=x)})
  SP          <- SpatialPolygons(Pols)
  SPDF        <- as(SP,"SpatialPolygonsDataFrame")

  #- Call the column to aggregate over
  if(is.null(control$clm)==TRUE){
    tacsat$count <- rep(1,nrow(tacsat))
    clm          <- "count"
  } else { clm <- control$clm }
  #- Do the aggregation
  tacsat$pol      <- poly_wkt
  if(is.null(control$FUN)==TRUE){
    agg             <- unlist(lapply(as.list(unique(tacsat$pol)),function(x){apply(t(tacsat[[clm]][which(tacsat$pol==x)]),1,sum,na.rm=TRUE)}))
  } else {
    agg             <- unlist(lapply(as.list(unique(tacsat$pol)),function(x){apply(t(tacsat[[clm]][which(tacsat$pol==x)]),1,control$FUN,na.rm=TRUE)}))
  }
  SPDF@data       <- data.frame(agg)
  colnames(SPDF@data) <- clm

  return(SPDF)
}

round2 <- function(x, n){
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
  }

polyDef <- function(lon, lat, gridx, gridy){
  # round to the nearest rectangle mid-point
  lon1 <- round2((lon - gridx/2)/gridx , 0) * gridx + gridx/2
  lat1 <- round2((lat - gridy/2)/gridy , 0) * gridy + gridy/2

  # create WKT sting
  out <- paste('Polygon(cbind(c('
    ,lon1 - gridx/2,','
    ,lon1 + gridx/2,','
    ,lon1 + gridx/2,','
    ,lon1 - gridx/2,','
    ,lon1 - gridx/2
    ,'),c('
    ,lat1 - gridy/2,','
    ,lat1 - gridy/2,','
    ,lat1 + gridy/2,','
    ,lat1 + gridy/2,','
    ,lat1 - gridy/2
    ,')))'
    ,sep='')

  return(out)
}


