#' Define a nested grid
#' 
#' Assign VMS points to a nested grid so that areas with a high density of
#' points will have small grid cells and areas with a low density will have
#' larger cells.
#' 
#' The alogrithm works as follows: A coarse starting grid is applied to the
#' positional data, the number of datapoints in each grid cell is counted and
#' if this number \code{>= n} then the cell is split in two.  Now the number of
#' datapoints in the smaller cells are counted again and any cells with
#' \code{>= n} will be split again, up to \code{maxstep} times.
#' 
#' This function allows data-rich areas to be plotted with a high spatial
#' resolution while allowing a lower spatial resolution for data-poor areas.
#' The nested grid also tends to reduce the amount of clustering within each
#' grid cell, which is important for estimating the area impacted by fishing
#' gear inside each cell. %% ~~ If necessary, more details than the description
#' above ~~
#' 
#' @param tacsat Tacsat dataframe
#' @param resx gridcell size in degrees in the longitude / x direction
#' @param resy gridcell size in degrees in the latitude / y direction
#' @param maxstep The maxiumum number of times the grid cells should be split
#' @param n If the number of points in a cell \code{>= n} then split the cell
#' @param control A list determining the output in the \code{data} slot, the
#' possible components are: \describe{ \item{list("clm")}{ The name of the
#' tacsat column that \code{FUN} should be applied to.  Also the name of the
#' only column in the dataframe that makes up the \code{data} slot.  The
#' default (\code{clm = NULL}) results \code{FUN} being applied a new column
#' called \code{"count"} which is simply \code{rep(1,nrow(tacsat))}.  }
#' \item{list("FUN")}{ The function to be applied to \code{tacsat[[clm]]}.  The
#' default (\code{FUN = NULL}) results in the function \code{sum}, so if
#' \code{clm = NULL} and \code{FUN = NULL} the result will be a count of the
#' number of datapoints.  } }
#' @return A SpatialPolygonsDataFrame, the value of the \code{data} slot will
#' be determined by the \code{control} settings.
#' @author Hans D. Gerritsen, Niels T. Hintzen
#' @seealso \code{\link{polyDef}}
#' @references Gerritsen, H. D., Minto, C. and Lordan, C. (2013) How much of
#' the seabed is impacted by mobile fishing gear? Absolute estimates from
#' Vessel Monitoring System (VMS) point data. ICES Journal of Marine Science
#' \bold{XX:XX}, doi: 10.1093/icesjms/fst017
#' @examples
#' 
#' data(tacsat)
#' tacsat            <- tacsat[sample(nrow(tacsat),2500),] # to speed it up a bit
#' tacsat            <- intervalTacsat(tacsat,level="vessel",fill.na=TRUE)
#' tacsat$INTV       <- ifelse(tacsat$INTV > 240, 240, tacsat$INTV)
#' tacsat$GEAR_WIDTH <- 0.024
#' tacsat$SWEPT_AREA <- tacsat$INTV / 60 * tacsat$SI_SP * tacsat$GEAR_WIDTH
#' SPDF              <- nestedGrid(tacsat, resx=1, resy=0.5, maxstep=10, n=20,
#'                       control=list(clm="SWEPT_AREA",FUN=sum))
#' SP                <- as(SPDF,"SpatialPolygons")
#' SP                <- surface(SP)
#' tempfun           <- function(x){lapply(x@Polygons,function(y){y@area})}
#' SPDF@data$surface <- unlist(lapply(SP@polygons,tempfun))
#' SPDF@data$SAratio <- SPDF@data$SWEPT_AREA / SPDF@data$surface
#' breaks            <- c(seq(0,quantile(SPDF@data$SAratio,0.95),length.out=9),2,5,10)
#' i                 <- cut(SPDF@data$SAratio,breaks=breaks)
#' SPDF@data$col     <- grey(seq(1, 0, length=length(breaks)))[i]
#' 
#' plot(NA, xlim = c(1, 5), ylim = c(51.5,55), xlab = 'Longitude',ylab = 'Latitude'
#'   ,asp=1/lonLatRatio(3,53))
#' plot(SP,col=SPDF@data$col,add=TRUE,border="lightgrey"); box()
#' points(tacsat$SI_LONG,tacsat$SI_LAT,cex=0.1,col=4)
#' 
#' @export nestedGrid
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
  SPDF        <- st_sf(agg,Pol)
  colnames(SPDF)[1] <- clm
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



#' Define polygon in text string format
#' 
#' Create a text string format polygon which can be called by other functions
#' 
#' The function stets up a grid - with the origin at (0,0) - and assigns a grid
#' cell to each of the points given by \code{lon} and \code{lat}.
#' 
#' @param lon Vector with longitude (or x) values
#' @param lat Vector with latitude (or y) values, should have the same length
#' as \code{lon}
#' @param gridx gridcell size in degrees in the longitude / x direction
#' @param gridy gridcell size in degrees in the latitude / y direction
#' @return A Well Known Text string for each value of \code{lon} and
#' \code{lat}.
#' @note Any points that lie exactly on the border between two grid cells will
#' be assigned to the grid cell above in the northern hemisphere, below in the
#' southern hemisphere, to the right in the eastern hemisphere and to the left
#' in the western hemisphere.
#' @author Hans D Gerritsen
#' @seealso \code{\link{nestedGrid}}
#' @examples
#' 
#' lon <- rnorm(25,3)
#' lat <- rnorm(25,53)
#' pols <- polyDef(lon, lat, gridx = 0.5, gridy= 0.5)
#' plot(lon,lat)
#' tempfun <- function(i){polygon(eval(parse(text=pols[i]))@coords)}
#' lapply(as.list(1:length(pols)),tempfun)
#' 
#' @export polyDef
polyDef <- function(lon, lat, gridx, gridy){
  # round to the nearest rectangle mid-point
  lon1 <- round2((lon - gridx/2)/gridx , 0) * gridx + gridx/2
  lat1 <- round2((lat - gridy/2)/gridy , 0) * gridy + gridy/2

  # create WKT sting
  out <- paste('st_polygon(list(cbind(c('
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
    ,'))))'
    ,sep='')

  return(out)
}


