
interpolation2Tacsat <- function(interpolation,tacsat,npoints=10,equalDist=TRUE){

# This function takes the list of tracks output by interpolateTacsat and converts them back to tacsat format.
# The npoints argument is the optional number of points between each 'real' position.
tacsat            <- sortTacsat(tacsat)
if(!"HL_ID" %in% colnames(tacsat)) tacsat$HL_ID <- 1:nrow(tacsat)
if(!"SI_DATIM" %in% colnames(tacsat)) tacsat$SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
if(equalDist){
  interpolationEQ <- equalDistance(interpolation,npoints)  #Divide points equally along interpolated track (default is 10).
} else {
  interpolationEQ <- lapply(interpolation,function(x){idx <- round(seq(2,nrow(x),length.out=npoints)); return(x[c(1,idx),])})
}


# Take off the tacsat row identifier
int       <- lapply(interpolationEQ,function(x){x[-1,]})
# Check if there are interpolations that were not succesful
lenint    <- lapply(int,length); idx <- which(unlist(lenint) == (npoints*2))
# Take off identifier again
int       <- lapply(interpolationEQ[idx],function(x){x[-1,]})
# Count number of interpolations
lenEQ     <- length(interpolationEQ[idx])
# get the interpolation
intmin2   <- lapply(int,function(x){x <- x[-c(1,nrow(x)),]; colnames(x)<- c("SI_LONG","SI_LATI"); return(x)})
# Get the identifiers
intidx    <- lapply(interpolationEQ[idx],function(x){x[1,]})
clnames   <- colnames(tacsat)
# Get the column names of the unique variables
b         <- lapply(intidx,function(x){cls <- clnames[tacsat[x[1],] %in% tacsat[x[2] ,]]; return(cls[!cls%in% c("SI_LONG","SI_LATI","SI_HE","SI_SP","SI_DATE","SI_TIME","SI_DATIM")])})
# make a matrix with new interpolated unique variables
bvals     <- lapply(as.list(1:lenEQ),function(x){matrix(unlist(tacsat[intidx[[x]][1],b[[x]]]),nrow=npoints-2,ncol=length(b[[x]]),byrow=TRUE,
                                                                 dimnames=list(round(seq(intidx[[x]][1],intidx[[x]][2],length.out=npoints-2),3),b[[x]]))})
# Create non-unique variables
SI_DATIMs <- data.frame(from=tacsat$SI_DATIM[do.call(rbind,intidx)[,1]],to=tacsat$SI_DATIM[do.call(rbind,intidx)[,2]])
SI_DATIMs <- lapply(as.list(1:lenEQ),function(x){seq(SI_DATIMs[x,1],SI_DATIMs[x,2],length.out=npoints)[2:(npoints-1)]})
SI_DATE   <- lapply(SI_DATIMs,function(x){format(x,format="%d/%m/%Y")})
timeNotation <- ifelse(length(unlist(strsplit(tacsat$SI_TIME[1],":")))>2,"secs","mins")
if(timeNotation=="secs") SI_TIME   <- lapply(SI_DATIMs,function(x){format(x,format="%H:%M:%S")})
if(timeNotation=="mins") SI_TIME   <- lapply(SI_DATIMs,function(x){format(x,format="%H:%M")})
SI_SPs    <- as.matrix(data.frame(from=tacsat$SI_SP[do.call(rbind,intidx)[,1]],to=tacsat$SI_SP[do.call(rbind,intidx)[,2]]))
SI_SP     <- mapply(seq,from=SI_SPs[,1],to=SI_SPs[,2],length.out=npoints-2)
SI_HE     <- lapply(as.list(1:lenEQ),function(x){y <- int[[x]]; return(bearing(y[2:(nrow(y)-1),1],y[2:(nrow(y)-1),2],y[3:nrow(y),1],y[3:nrow(y),2]))})
HL_ID     <- lapply(as.list(1:lenEQ),function(x){rep(tacsat$HL_ID[intidx[[x]][1]],npoints-2)})

# Combine all
ret       <- lapply(as.list(1:lenEQ),function(x){data.frame(bvals[[x]],
                                                            intmin2[[x]],
                                                            SI_DATIM=SI_DATIMs[[x]],
                                                            SI_DATE=SI_DATE[[x]],
                                                            SI_TIME=SI_TIME[[x]],
                                                            SI_SP=SI_SP[,x],
                                                            SI_HE=SI_HE[[x]],
                                                            HL_ID=HL_ID[[x]],stringsAsFactors=F)})
# Remove duplicate records
idx       <- which(unlist(lapply(ret,function(x){length(grep(".1",colnames(x)))>0}))==TRUE)
nidx      <- which(unlist(lapply(ret,function(x){length(grep(".1",colnames(x)))>0}))==FALSE)
retmin    <- lapply(ret[idx],function(x){x[,-grep(".1",colnames(x))]})
ret       <- c(retmin,ret[nidx])

# Order columns
ret       <- lapply(ret,function(x){addclnames <- clnames[which(!clnames %in% colnames(x))];
                                    x[,addclnames] <- NA;
                                    x[,clnames]})
# Combine all interpolations with original file and return
interpolationTot <- do.call(rbind,ret)
interpolationTot <- formatTacsat(interpolationTot)
tacsatInt <- rbindTacsat(tacsat,interpolationTot)
tacsatInt <- sortTacsat(tacsatInt)

return(tacsatInt)}

                                                                #467