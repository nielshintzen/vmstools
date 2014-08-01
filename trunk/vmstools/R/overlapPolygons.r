# overlapPolygons  to calculate the intersection area of 2 or more polygons

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
