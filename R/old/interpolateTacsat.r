`interpolateTacsat` <-
function(tacsat                          #VMS datapoints
                              ,interval=120             #Specify in minutes, NULL means use all points
                              ,margin=12                #Specify the margin in minutes that the interval might deviate in a search for the next point
                              ,res=100                  #Resolution of interpolation method (default = 100)
                              ,method="cHs"             #Specify the method to be used: Straight line (SL) of cubic Hermite spline (cHs)
                              ,params=list(fm=0.5,distscale=20,sigline=0.2,st=c(2,6))  #Specify the three parameters: fm, distscale, sigline, speedthreshold
                              ,headingAdjustment=0
                              ,fast=FALSE){

if(!"SI_DATIM" %in% colnames(tacsat)) tacsat$SI_DATIM     <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
                              
  #Start interpolating the data
if(!method %in% c("cHs","SL"))  stop("method selected that does not exist")

#-------------------------------------------------------------------------------
#Fast method or not
#-------------------------------------------------------------------------------
if(fast){
  #Interpolation only by vessel, so split tacsat up
  tacsat$ID <- 1:nrow(tacsat)
  splitTa   <- split(tacsat,tacsat$VE_REF)
  spltTaCon <- lapply(splitTa,function(spltx){
                  #Calculate time different between every record
                  dftimex <- outer(spltx$SI_DATIM,spltx$SI_DATIM,difftime,units="mins")
                  iStep   <- 1
                  connect <- list()
                  counter <- 1
                  #Loop over all possible combinations and store if a connection can be made
                  while(iStep <= nrow(spltx)){
                    endp <- which(dftimex[,iStep] >= (interval - margin) & dftimex[,iStep] <= (interval + margin))
                    if(length(endp)>0){
                      if(length(endp)>1) endp <- endp[which.min(abs(interval - dftimex[endp,iStep]))][1]
                      connect[[counter]]    <- c(iStep,endp)
                      counter               <- counter + 1
                      iStep                 <- endp
                    } else { iStep          <- iStep + 1}
                  }
                  #Return matrix of conenctions
                  return(do.call(rbind,connect))})

  if(method=="cHs") returnInterpolations <- unlist(lapply(as.list(names(unlist(lapply(spltTaCon,nrow)))),function(y){
                                              return(interCubicHermiteSpline(spltx=splitTa[[y]],spltCon=spltTaCon[[y]],res,params,headingAdjustment))}),recursive=FALSE)
  if(method=="SL")  returnInterpolations <- unlist(lapply(as.list(names(unlist(lapply(spltTaCon,nrow)))),function(y){
                                              return(interStraightLine(splitTa[[y]],spltTaCon[[y]],res))}),recursive=FALSE)

} else {
  

    #Initiate returning result object
  returnInterpolations <- list()

    #Start iterating over succeeding points
  for(iStep in 1:(dim(tacsat)[1]-1)){
    if(iStep == 1){
      iSuccess    <- 0
      endDataSet  <- 0
      startVMS    <- 1
      ship        <- tacsat$VE_REF[startVMS]
    } else {
        if(is.na(endVMS)==TRUE) endVMS <- startVMS + 1
        startVMS <- endVMS
        #-Check if the end of the dataset is reached
        if(endDataSet == 1 & rev(unique(tacsat$VE_REF))[1] != ship){
          startVMS  <- which(tacsat$VE_REF == unique(tacsat$VE_REF)[which(unique(tacsat$VE_REF)==ship)+1])[1]
          ship      <- tacsat$VE_REF[startVMS]
          endDataSet<- 0
        }
        if(endDataSet == 1 & rev(unique(tacsat$VE_REF))[1] == ship) endDataSet <- 2 #Final end of dataset
      }

    #if end of dataset is not reached, try to find succeeding point
    if(endDataSet != 2){
      result      <- findEndTacsat(tacsat,startVMS,interval,margin)
      endVMS      <- result[1]
      endDataSet  <- result[2]
      if(is.na(endVMS)==TRUE) int <- 0  #No interpolation possible
      if(is.na(endVMS)==FALSE) int <- 1  #Interpolation possible

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

        if (is.na(tacsat[startVMS,"SI_HE"])=="TRUE") tacsat[startVMS,"SI_HE"] <- 0
        if (is.na(tacsat[endVMS,  "SI_HE"])=="TRUE") tacsat[endVMS,  "SI_HE"] <- 0

          #Heading at begin point in degrees
        Hx0 <- sin(tacsat[startVMS,"SI_HE"]/(180/pi))
        Hy0 <- cos(tacsat[startVMS,"SI_HE"]/(180/pi))
          #Heading at end point in degrees
        Hx1 <- sin(tacsat[endVMS-headingAdjustment,"SI_HE"]/(180/pi))
        Hy1 <- cos(tacsat[endVMS-headingAdjustment,"SI_HE"]/(180/pi))

        Mx0 <- tacsat[startVMS, "SI_LONG"]
        Mx1 <- tacsat[endVMS,   "SI_LONG"]
        My0 <- tacsat[startVMS, "SI_LATI"]
        My1 <- tacsat[endVMS,   "SI_LATI"]

          #Corrected for longitude lattitude effect
        Hx0 <- Hx0 * params$fm * tacsat[startVMS,"SI_SP"] /((params$st[2]-params$st[1])/2+params$st[1])
        Hx1 <- Hx1 * params$fm * tacsat[endVMS,"SI_SP"]   /((params$st[2]-params$st[1])/2+params$st[1])
        Hy0 <- Hy0 * params$fm * lonLatRatio(tacsat[c(startVMS,endVMS),"SI_LONG"],tacsat[c(startVMS,endVMS),"SI_LATI"])[1] * tacsat[startVMS,"SI_SP"]/((params$st[2]-params$st[1])/2+params$st[1])
        Hy1 <- Hy1 * params$fm * lonLatRatio(tacsat[c(startVMS,endVMS),"SI_LONG"],tacsat[c(startVMS,endVMS),"SI_LATI"])[2] * tacsat[endVMS,"SI_SP"]/((params$st[2]-params$st[1])  /2+params$st[1])

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
        fx <- seq(tacsat$SI_LONG[startVMS],tacsat$SI_LONG[endVMS],length.out=res)
        fy <- seq(tacsat$SI_LATI[startVMS],tacsat$SI_LATI[endVMS],length.out=res)

          #Add one to list of successful interpolations
        iSuccess <- iSuccess + 1
        returnInterpolations[[iSuccess]] <- matrix(rbind(c(startVMS,endVMS),cbind(fx,fy)),ncol=2,dimnames=list(c("startendVMS",seq(1,res,1)),c("x","y")))
      }
    }
  }
}
  
return(returnInterpolations)}

