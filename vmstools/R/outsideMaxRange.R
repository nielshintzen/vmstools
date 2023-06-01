#' compute fraction of Confidence Interval that is located outside a maximum
#' range
#' 
#' The calculation of the Confidence Interval surrounding an interpolation
#' depends on two parameters: sigline & distscale. These use of these
#' parameters could result in extremely wide or extremely small CI's. To check
#' which proportion is located inside and outside the maximum range (as defined
#' by an ellipse), this function calculates this proportion, as well as the
#' maximum value representing the starting and end point CI values.
#' 
#' 
#' @param int interpolation, as data.frame from 'interpolateTacsat'
#' @param tacint tacsat records (two rows) corresponding with the interpolation
#' @param params list of parameters used to perform interpolation
#' @param grid object of class 'GridTopology' specifying the grid dimensions
#' @return Returnes a list with three objects: 1) Fraction of CI located
#' outside maximum range 2) Franction of CI located inside maximum range 3)
#' Maximum value of CI (top, at VMS points)
#' @author Niels T. Hintzen
#' @seealso \code{\link{createGrid}}, \code{\link{calculateCI}},
#' \code{\link{point.in.polygon}}
#' @references Hintzen et al. 2010 Improved estimation of trawling tracks using
#' cubic Hermite spline interpolation of position registration data, EU lot 2
#' project
#' @examples
#' 
#' data(tacsat)
#' 
#'   #Sort the Tacsat data
#' tacsat     <- sortTacsat(tacsat)
#' tacsat     <- tacsat[1:1000,]
#' 
#'   #Filter the Tacsat data
#' tacsat          <- filterTacsat(tacsat,c(2,6),hd=NULL,remDup=TRUE)
#' 
#'   #Interpolate the VMS data
#' interpolation <- interpolateTacsat(tacsat,interval=120,margin=10,
#'                     res=100,method="cHs",params=list(fm=0.5,distscale=20,
#'                     sigline=0.2,st=c(2,6)),headingAdjustment=0)
#' 
#'   #Create the final grid where all interpolations should fit on
#' xrange        <- c(2,3); yrange <- c(51,52)
#' grid          <- createGrid(xrange,yrange,resx=0.01,resy=0.005)
#' 
#' res           <- outsideMaxRange(interpolation[[4]],
#'                                  tacsat[interpolation[[4]][1,],],
#'                                  params=list(fm=0.25,distscale=3.1,
#'                                              sigline=0.4,st=c(2,6)),
#'                                  grid=grid)
#' 
#' 
#' @export outsideMaxRange
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
