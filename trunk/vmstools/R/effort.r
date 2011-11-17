effort <- function(x,level="trip",unit="hours",weight=c(0.5,0.5),fill.na=F){

              if("SI_FT" %in% colnames(x)) x$FT_REF <- x$SI_FT
              #-Add if necessary a datim column
              if(!"SI_DATIM" %in% colnames(x)){
                if(all(c("VE_FLT","VE_KW") %in% colnames(x)))    x$SI_DATIM <- as.POSIXct(paste(x$FT_DDAT,  x$FT_DTIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
                if(all(c("SI_LATI","SI_LONG") %in% colnames(x))) x$SI_DATIM <- as.POSIXct(paste(x$SI_DATE,  x$SI_TIME,    sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
              }

              if(level != "trip") stop("No other level than trip possible")
              if(!length(grep("FT_REF",colnames(x)))>0) stop("No FT_REF defined")

              #- Remove trips without trip number
              x <- subset(x,FT_REF != 0)

              #- Order the data
              x$ID  <- paste(x$VE_REF,x$FT_REF,sep="_")
              x     <- orderBy(~VE_REF+SI_DATIM+FT_REF,data=x)

              if(all(c("SI_LATI","SI_LONG") %in% colnames(x))){
                x$LE_EFF_VMS  <- intervalTacsat(x,level="trip",weight=weight,fill.na=fill.na)
                if(!unit %in% c("secs","mins","hours","days","weeks")) stop("Unit must be in 'secs,mins,hours,days or weeks'")
                if(unit == "secs")    x$LE_EFF_VMS  <- x$LE_EFF_VMS * 60
                if(unit == "hours")   x$LE_EFF_VMS  <- x$LE_EFF_VMS/60
                if(unit == "days")    x$LE_EFF_VMS  <- x$LE_EFF_VMS / 60 / 24
                if(unit == "weeks")   x$LE_EFF_VMS  <- x$LE_EFF_VMS / 60 / 24 / 7
              }

              if(all(c("VE_FLT","VE_KW") %in% colnames(x))){
                x$SI_DATIM2   <- as.POSIXct(paste(x$FT_LDAT,  x$FT_LTIME,    sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
                x$LE_EFF_LOG  <- an(difftime(x$SI_DATIM2,x$SI_DATIM,units=unit))
              }

        return(x)}