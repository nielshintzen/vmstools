`interpolateVMS` <-
function(tacsat                          #VMS datapoints
                              ,interval=120             #Specify in minutes, NULL means use all points
                              ,margin=12                #Specify the margin in minutes that the interval might deviate in a search for the next point
                              ,res=100                  #Resolution of interpolation method (default = 100)
                              ,method="cHs"             #Specify the method to be used: Straight line (SL) of cubic Hermite spline (cHs)
                              ,params=list(fm=0.5,distscale=20,sigline=0.2,st=c(2,6))  #Specify the three parameters: fm, distscale, sigline, speedthreshold
                              ,headingAdjustment=0
                              ){

VMS. <- tacsat
colnames(VMS.) <- c("ship","lat","lon","date","time","speed","heading")
VMS.$datim     <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S")
                              
  #Start interpolating the data
if(!method %in% c("cHs","SL"))  stop("method selected that does not exist")

  #Initiate returning result object
returnInterpolations <- list()

  #Start iterating over succeeding points
for(iStep in 1:(dim(VMS.)[1]-1)){
  if(iStep == 1){
    iSuccess    <- 0
    endDataSet  <- 0
    startVMS    <- 1
    ship        <- VMS.$ship[startVMS]
  } else {
      if(is.na(endVMS)==T) endVMS <- startVMS + 1
      startVMS <- endVMS
      if(endDataSet == 1 & rev(unique(VMS.$ship))[1] != ship){
        startVMS  <- which(VMS.$ship == unique(VMS.$ship)[which(unique(VMS.$ship)==ship)+1])[1]
        ship      <- VMS.$ship[startVMS]
        endDataSet<- 0
      }
      if(endDataSet == 1 & rev(unique(VMS.$ship))[1] == ship) endDataSet <- 2 #Final end of dataset
    }


  if(endDataSet != 2){
    result      <- findEndTacsat(VMS.,startVMS,interval,margin)
    endVMS      <- result[1]
    endDataSet  <- result[2]
    if(is.na(endVMS)==T) int <- 0  #No interpolation possible
    if(is.na(endVMS)==F) int <- 1  #Interpolation possible
  
      #Interpolate according to the Cubic Hermite Spline method
    if(method == "cHs" & int == 1){

        #Define the cHs formula
      F00 <- numeric()
      F10 <- numeric()
      F01 <- numeric()
      F11 <- numeric()
      i   <- 0
      t   <- seq(0,1,length.out=res)
      F00 <- 2*t^3 -3*t^2 + 1
      F10 <- t^3-2*t^2+t
      F01 <- -2*t^3+3*t^2
      F11 <- t^3-t^2          

      if (is.na(VMS.[startVMS,"heading"])=="TRUE") VMS.[startVMS,"heading"] <- 0
      if (is.na(VMS.[endVMS,  "heading"])=="TRUE") VMS.[endVMS,  "heading"] <- 0
      
        #Heading at begin point in degrees
      Hx0 <- sin(VMS.[startVMS,"heading"]/(180/pi))
      Hy0 <- cos(VMS.[startVMS,"heading"]/(180/pi))
        #Heading at end point in degrees
      Hx1 <- sin(VMS.[endVMS-headingAdjustment,"heading"]/(180/pi))
      Hy1 <- cos(VMS.[endVMS-headingAdjustment,"heading"]/(180/pi))
      
      Mx0 <- VMS.[startVMS, "declon"]
      Mx1 <- VMS.[endVMS,   "declon"]
      My0 <- VMS.[startVMS, "declat"]
      My1 <- VMS.[endVMS,   "declat"]

        #Corrected for longitude lattitude effect
      Hx0 <- Hx0 * params$fm * VMS.[startVMS,"speed"] /((params$st[2]-params$st[1])/2+params$st[1])
      Hx1 <- Hx1 * params$fm * VMS.[endVMS,"speed"]   /((params$st[2]-params$st[1])/2+params$st[1])
      Hy0 <- Hy0 * params$fm * lonLatRatio(VMS.[c(startVMS,endVMS),"declon"],VMS.[c(startVMS,endVMS),"declat"])[1] * VMS.[startVMS,"speed"]/((params$st[2]-params$st[1])/2+params$st[1])
      Hy1 <- Hy1 * params$fm * lonLatRatio(VMS.[c(startVMS,endVMS),"declon"],VMS.[c(startVMS,endVMS),"declat"])[2] * VMS.[endVMS,"speed"]/((params$st[2]-params$st[1])  /2+params$st[1])

        #Finalizing the interpolation based on cHs
      fx <- numeric()
      fy <- numeric()
      fx <- F00*Mx0+F10*Hx0+F01*Mx1+F11*Hx1
      fy <- F00*My0+F10*Hy0+F01*My1+F11*Hy1
      
        #Add one to list of successful interpolations
      iSuccess <- iSuccess + 1 
      returnInterpolations[[iSuccess]] <- matrix(rbind(c(startVMS,endVMS),cbind(fx,fy)),ncol=2,dimnames=list(c("startendVMS",seq(1,res,1)),c("x","y")))
    }
    
      #Interpolate according to a straight line
    if(method == "SL" & int == 1){
      fx <- seq(VMS.$declon[startVMS],VMS.$declon[endVMS],length.out=res)
      fy <- seq(VMS.$declat[startVMS],VMS.$declat[endVMS],length.out=res)
      
        #Add one to list of successful interpolations
      iSuccess <- iSuccess + 1
      returnInterpolations[[iSuccess]] <- matrix(rbind(c(startVMS,endVMS),cbind(fx,fy)),ncol=2,dimnames=list(c("startendVMS",seq(1,res,1)),c("x","y")))
    }
  }
}
return(returnInterpolations)}

