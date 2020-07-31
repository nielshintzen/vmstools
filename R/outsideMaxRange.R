outsideMaxRange <- function(int
                               ,tacint
                               ,params
                               ,grid){
    if (!"SI_DATIM" %in% colnames(tacint))
        tacint$SI_DATIMIM <- as.POSIXct(paste(tacint$SI_DATE, tacint$SI_TIME,
            sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M")
                               
    CI      <- calculateCI(int,
                               tacint,
                               params,
                               grid,
                               plot=FALSE)

    mxr     <-  maxRangeCI(x=c(int[-1,1][1],rev(int[-1,1])[1]),
                           y  =c(int[-1,2][1],rev(int[-1,2])[1]),
                           Time=c(difftime(tacint$SI_DATIM[2],tacint$SI_DATIM[1],units="mins")),
                           speed=pmax(tacint$SI_SP,rep(distanceInterpolation(list(int)) / 1.852 /
                                      c(difftime(tacint$SI_DATIM[2],tacint$SI_DATIM[1],units="hours")),2)))

    coords  <- coordinates(CI)
    propCI  <- point.in.polygon(coords[,1],coords[,2],mxr[[1]][,1],mxr[[1]][,2])
    insideR <- sum(CI@data$data[which(propCI == 1)],na.rm=TRUE) / sum(CI@data$data,na.rm=TRUE) #Sum of total CI values inside the maximum range, should ideally be all the grid cells with values
    outsideR<- sum(CI@data$data[which(propCI == 0)],na.rm=TRUE) / sum(CI@data$data,na.rm=TRUE)#Sum of total CI values outside the maximum range, should ideally be 0
    maxR    <- max(CI@data$data[which(propCI == 1)],na.rm=TRUE) #Top of the CI, should ideally equal to 1
return(list(insideR,outsideR,maxR))}
