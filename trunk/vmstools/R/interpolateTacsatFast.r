`interpolateTacsatFast` <-
function(tacsat                                         #VMS datapoints
                              ,interval=120             #Specify in minutes, NULL means use all points
                              ,margin=12                #Specify the margin in minutes that the interval might deviate in a search for the next point
                              ,res=100                  #Resolution of interpolation method (default = 100)
                              ,method="cHs"             #Specify the method to be used: Straight line (SL) of cubic Hermite spline (cHs)
                              ,params=list(fm=0.5,distscale=20,sigline=0.2,st=c(2,6))  #Specify the three parameters: fm, distscale, sigline, speedthreshold
                              ,headingAdjustment=0
                              ){

if(!"SI_DATIM" %in% colnames(tacsat)) tacsat$SI_DATIM     <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
tacsat$ID <- 1:nrow(tacsat) #need ID for

  #Start interpolating the data
if(!method %in% c("cHs","SL"))  stop("method selected that does not exist")

  #Interpolation only by vessel, so split tacsat up
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
                
if(method=="cHs") spltTaInt <- unlist(lapply(as.list(names(unlist(lapply(spltTaCon,nrow)))),function(y){
                                  return(interCubicHermiteSpline(spltx=splitTa[[y]],spltCon=spltTaCon[[y]],res,params,headingAdjustment))}),recursive=F)
if(method=="SL")  spltTaInt <- unlist(lapply(as.list(names(unlist(lapply(spltTaCon,nrow)))),function(y){
                                  return(interStraightLine(splitTa[[y]],spltTaCon[[y]],res))}),recursive=F)

return(spltTaInt)}
                
