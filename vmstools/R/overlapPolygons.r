# overlapPolygons  to calculate the intersection area of 2 or more polygons



#' Calculate the surface area that polygons have in common
#' 
#' Calculate the surface area that 2 or more polygons have in common
#' 
#' 
#' @param pol1 Polygon number 1. Can be of class 'data.frame' with first column
#' Longitude and second column Latitude. Can be of class 'PolySet' and can be
#' of class 'SpatialPolygons'.
#' @param pol2 Polygon number 2. Can be of class 'data.frame' with first column
#' Longitude and second column Latitude. Can be of class 'PolySet' and can be
#' of class 'SpatialPolygons'.
#' @param projection Optional projection attribute to add (often "LL" for
#' Longitude-Latitude).
#' @param zone Optional zone attribute to add.
#' @return Returns a data.frame with overlap in km2 (or non-projected units)
#' with PID referring to the combination of the two polygon sets.
#' @author Katell Hamon, Niels Hintzen
#' @seealso \code{\link{lonLat2SpatialPolygons}}, \code{\link{as.PolySet}},
#' \code{\link{surface}}
#' @examples
#' 
#' 
#' 
#' #- Test with data.frame polygons
#' pol1 <- data.frame(cbind(c(2,3,3,2),c(54,54,55,55)))
#' pol2 <- data.frame(cbind(c(2,3,3,2),c(55,55,54.5,54.5)))
#' 
#' overlapPolygons(pol1,pol2)
#' 
#' #- Test with SpatialPolygons
#' pol1 <- lonLat2SpatialPolygons(SI_LONG=c(2,3,3,2),SI_LATI=c(54,54,55,55))
#' pol2 <- lonLat2SpatialPolygons(SI_LONG=c(2,3,3,2),SI_LATI=c(54.5,54.5,55,55))
#' 
#' overlapPolygons(pol1,pol2)
#' surface(pol1)@polygons[[1]]@Polygons[[1]]@area
#' 
#' #- Test with PolySet polygons
#' pol1 <- as.PolySet(data.frame(PID=rep(1,4),POS=1:4,X=c(2,3,3,2),Y=c(54,54,55,55)))
#' pol2 <- as.PolySet(data.frame(PID=rep(1,4),POS=1:4,X=c(2,3,3,2),Y=c(54.5,54.5,55,55)))
#' 
#' overlapPolygons(pol1,pol2)
#' 
#' #- Test with multiple polygons
#' data(tacsat)
#' pols1 <- cbind(s1=rep(2,length(seq(49,63,2))),s2=c(seq(49,63,2)))
#' pols2 <- cbind(s1=tacsat$SI_LONG[seq(2,nrow(tacsat),length.out=5)],s2=tacsat$SI_LATI[seq(2,nrow(tacsat),length.out=5)])
#' resx   <- 1; resy <- 0.5
#' sPols1                       <- lonLat2SpatialPolygons(lst=lapply(as.list(1:nrow(pols1)),
#'                                   function(x){data.frame(SI_LONG=c(pols1[x,"s1"]-resx/2,rep(pols1[x,"s1"]+resx/2,2),pols1[x,"s1"]-resx/2),
#'                                                          SI_LATI=c(rep(pols1[x,"s2"]-resy/2,2),rep(pols1[x,"s2"]+resy/2,2)))}))
#' sPols2                       <- lonLat2SpatialPolygons(lst=lapply(as.list(1:nrow(pols2)),
#'                                   function(x){data.frame(SI_LONG=c(pols2[x,"s1"]-resx/2,rep(pols2[x,"s1"]+resx/2,2),pols2[x,"s1"]-resx/2),
#'                                                          SI_LATI=c(rep(pols2[x,"s2"]-resy/2,2),rep(pols2[x,"s2"]+resy/2,2)))}))
#' 
#' overlapPolygons(sPols1,sPols2)
#' 
#' 
#' @export overlapPolygons
overlapPolygons <- function(pol1=NULL,pol2=NULL,projection="LL",zone=NULL){

  #- Class = Polyset
  if(class(pol1)[1]=="PolySet")
    pol1 <- as.PolySet(pol1,projection=projection,zone=zone)
  if(class(pol2)[1]=="PolySet")
    pol2 <- as.PolySet(pol2,projection=projection,zone=zone)

  #- Class = data.frame
  if(class(pol1)[1] == "data.frame"){
    if(nrow(pol1)>ncol(pol1))
      pol1 <- as.PolySet(data.frame(PID=1,POS=1:nrow(pol1),X=pol1[,1],Y=pol1[,2]),projection=projection,zone=zone)
    if(nrow(pol1)<ncol(pol1))
      pol1 <- as.PolySet(data.frame(PID=1,POS=1:nrow(pol1),X=pol1[1,],Y=pol1[2,]),projection=projection,zone=zone)
   }
   if(class(pol2)[1] == "data.frame"){
    if(nrow(pol2)>ncol(pol2))
      pol2 <- as.PolySet(data.frame(PID=1,POS=1:nrow(pol2),X=pol2[,1],Y=pol2[,2]),projection=projection,zone=zone)
    if(nrow(pol2)<ncol(pol2))
      pol2 <- as.PolySet(data.frame(PID=1,POS=1:nrow(pol2),X=pol2[1,],Y=pol2[2,]),projection=projection,zone=zone)
   }
   
   #- Class = SpatialPolygons
   if(class(pol1)[1]=="SpatialPolygons"){
     pols1 <- numeric()
     counter       <- 0
     for(iPol1 in 1:length(pol1@polygons)){
       for(iPol2 in 1:length(pol1@polygons[[iPol1]]@Polygons)){
         counter   <- counter + 1
         sourcePoly <- data.frame(cbind(1,1:nrow(pol1@polygons[[iPol1]]@Polygons[[iPol2]]@coords),
                                  pol1@polygons[[iPol1]]@Polygons[[iPol2]]@coords[,1],pol1@polygons[[iPol1]]@Polygons[[iPol2]]@coords[,2]))
         rownames(sourcePoly)<-1:nrow(sourcePoly)
         colnames(sourcePoly)<-c("PID","POS","X","Y")
         sourcePoly$PID[]    <-counter

         pols1       <- rbind(pols1,sourcePoly)
       }
     }

     pol1   <- as.PolySet(data.frame(pols1),projection=projection,zone=zone)
   }
   if(class(pol2)[1]=="SpatialPolygons"){
     pols2 <- numeric()
     counter       <- 0
     for(iPol1 in 1:length(pol2@polygons)){
       for(iPol2 in 1:length(pol2@polygons[[iPol1]]@Polygons)){
         counter   <- counter + 1
         sourcePoly <- data.frame(cbind(1,1:nrow(pol2@polygons[[iPol1]]@Polygons[[iPol2]]@coords),
                                  pol2@polygons[[iPol1]]@Polygons[[iPol2]]@coords[,1],pol2@polygons[[iPol1]]@Polygons[[iPol2]]@coords[,2]))
         rownames(sourcePoly)<-1:nrow(sourcePoly)
         colnames(sourcePoly)<-c("PID","POS","X","Y")
         sourcePoly$PID[]    <-counter

         pols2       <- rbind(pols2,sourcePoly)
       }
     }
    pol2   <- as.PolySet(data.frame(pols2),projection=projection,zone=zone)
   }


  jpoly  <- try(joinPolys(pol1,pol2))
  if (!is.null(jpoly) & class(jpoly)!= "try-error"){
    ar <- calcArea(jpoly)
    jarea  <- aggregate(ar['area'],by=list(PID=ar$PID),FUN=sum)
  }else{
    jarea <- data.frame(PID=unique(pol1$PID),area=0)
  }
  return(jarea)
}
