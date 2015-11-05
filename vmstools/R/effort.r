effort <- function(x,by="FT_REF",unit="hours",weight=c(0.5,0.5),fill.na=FALSE){
              dattype <- ifelse(all(c("SI_LATI","SI_LONG") %in% colnames(x)),"tacsat","eflalo")

              #-Add if necessary a datim column
              if(!"SI_DATIM" %in% colnames(x)){
                if(dattype=="eflalo") x$SI_DATIM <- as.POSIXct(paste(x$FT_DDAT,  x$FT_DTIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
                if(dattype=="tacsat") x$SI_DATIM <- as.POSIXct(paste(x$SI_DATE,  x$SI_TIME,    sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
              }

              if(!length(grep("FT_REF",colnames(x)))>0 & "FT_REF" %in% by)   stop("No trip reference (FT_REF) defined. For tacsat, use mergeEflalo2Tacsat() first")
              if(!length(grep("LE_GEAR",colnames(x)))>0 & "LE_GEAR" %in% by) stop("No gear (LE_GEAR) defined. For tacsat, use mergeEflalo2Tacsat() first")
              if("byRow"%in% by & length(by)>1) stop("If 'byRow' is selected, no other attributed can be added")

              #- Remove trips without trip number
              x <- subset(x,FT_REF != 0)

              #- Order the data
              x$ID  <- paste(x$VE_REF,x$FT_REF,sep="_")
              x     <- orderBy(~VE_REF+SI_DATIM+FT_REF,data=x)

              #- Add date notation
              if("SI_DAY" %in% by){
                if(dattype == "eflalo") x$SI_DAY <- yday(as.POSIXct(x$LE_CDAT,format="%d/%m/%Y"))
                if(dattype == "tacsat") x$SI_DAY <- yday(x$SI_DATIM)
              }
              if("SI_YEAR" %in% by){
                if(dattype == "eflalo") x$SI_YEAR<- year(as.POSIXct(x$LE_CDAT,format="%d/%m/%Y"))
                if(dattype == "tacsat") x$SI_YEAR<- year(x$SI_DATIM)
              }
              if("SI_MONTH" %in% by){
                if(dattype == "eflalo") x$SI_MONTH<- month(as.POSIXct(x$LE_CDAT,format="%d/%m/%Y"))
                if(dattype == "tacsat") x$SI_MONTH<- month(x$SI_DATIM)
              }
              if("SI_WEEK" %in% by){
                if(dattype == "eflalo") x$SI_WEEK<- week(as.POSIXct(x$LE_CDAT,format="%d/%m/%Y"))
                if(dattype == "tacsat") x$SI_WEEK<- week(x$SI_DATIM)
              }
              if("SI_QUARTER" %in% by){
                if(dattype == "eflalo") x$SI_QUARTER<- quarter(as.POSIXct(x$LE_CDAT,format="%d/%m/%Y"))
                if(dattype == "tacsat") x$SI_QUARTER<- quarter(x$SI_DATIM)
              }

              if("LE_RECT" %in% by)
                if(dattype == "tacsat") x$LE_RECT <- ICESrectangle(x)
              if("LE_AREA" %in% by){
                data(ICESareas)
                if(dattype == "tacsat"){
                  x$LE_AREA <- ICESarea(x,ICESareas)
                  x$LE_AREA[which(is.na(x$LE_AREA)==T)] <- "OTHER"
                }
                if(dattype == "eflalo"){
                  eflonlat        <- ICESrectangle2LonLat(x$LE_RECT)
                  x$LE_AREA       <- ICESarea(eflonlat,ICESareas)
                  x$LE_AREA[which(is.na(x$LE_AREA)==T)] <- "OTHER"
                }
              }
              #-----------------------------------------------------------------------------
              #- Calculate effort
              #-----------------------------------------------------------------------------
              if(dattype == "eflalo"){
                  if(!"FT_DDATIM" %in% colnames(x))
                    x$FT_DDATIM <- as.POSIXct(paste(x$FT_DDAT,x$FT_DTIME),format="%d/%m/%Y %H:%M",tz = "GMT")
                  if(!"FT_LDATIM" %in% colnames(x))
                    x$FT_LDATIM <- as.POSIXct(paste(x$FT_LDAT,x$FT_LTIME),format="%d/%m/%Y %H:%M",tz = "GMT")

                x$LE_CDATIM <- as.POSIXct(x$LE_CDAT,format="%d/%m/%Y",tz="GMT")
                x$INTV      <- c(difftime(x$FT_LDATIM,x$FT_DDATIM,units="mins"))
                x$FT_DURDAY <- ifelse(c(difftime(as.Date(x$FT_LDATIM),as.Date(x$FT_DDATIM),units="hours") == 0),
                                      c(difftime(x$FT_LDATIM,x$FT_DDATIM,units="mins")),
                                      ifelse(c(difftime(as.Date(x$FT_DDATIM),as.Date(x$LE_CDATIM),units="hours")==0),
                                        c(difftime(x$LE_CDATIM+(60*60*24),x$FT_DDATIM,units="mins")),
                                        ifelse(c(difftime(as.Date(x$FT_LDATIM),as.Date(x$LE_CDATIM),units="hours")==0),
                                          c(difftime(x$FT_LDATIM,x$LE_CDATIM,units="mins")),
                                          1440)))
                # Here there is still a problem because INTVDAY is calculated for catch days only, so you miss some effort of a whole trip
                x$dummy     <- 1
                x           <- merge(x,aggregate(x$dummy,by=list(x$FT_REF,x$LE_CDATIM),FUN=sum,na.rm=T),by.x=c("FT_REF","LE_CDATIM"),by.y=c("Group.1","Group.2"),all.x=T)
                colnames(x)[length(colnames(x))] <- "NR_FT_REF"
                if("day" %in% by){
                  x$INTVDAY   <- x$FT_DURDAY / x$NR_FT_REF
                } else {
                  x$INTVDAY   <- x$INTV / x$NR_FT_REF
                  }
                x$INTV      <- x$INTVDAY
                x           <- x[,-grep("INTVDAY",colnames(x))]
                x           <- x[,-grep("dummy",colnames(x))]
                x           <- x[,-grep("FT_DURDAY",colnames(x))]
                x           <- x[,-grep("NR_FT_REF",colnames(x))]
              }
              if(dattype == "tacsat"){
                  x$INTV  <- intervalTacsat(x,level="vessel",weight=weight,fill.na=fill.na)$INTV
                #- Overwrite INTV if also trip in level
                if("FT_REF" %in% by | "FT_REF" %in% colnames(x))
                  x$INTV  <- intervalTacsat(x,level="trip",weight=weight,fill.na=fill.na)$INTV
               }
               if(!"byRow" %in% by){
                 if(length(by)==1) res        <- aggregate(x$INTV,by=list(x[,by]),FUN=sum,na.rm=T)
                 if(length(by)!=1) res        <- aggregate(x$INTV,by=as.list(x[,by]),FUN=sum,na.rm=T)
                 colnames(res)                   <- c(by,"EFFORT")
               }

               if("byRow" %in% by){
                 x$EFFORT   <- x$INTV
                 if(unit == "secs")    x$EFFORT  <- x$EFFORT * 60
                 if(unit == "hours")   x$EFFORT  <- x$EFFORT/60
                 if(unit == "days")    x$EFFORT  <- x$EFFORT / 60 / 24
                 if(unit == "weeks")   x$EFFORT  <- x$EFFORT / 60 / 24 / 7
               }
               if(!"byRow" %in% by){
                 x                               <- res
                 if(!unit %in% c("secs","mins","hours","days","weeks")) stop("Unit must be in 'secs,mins,hours,days or weeks'")
                 if(unit == "secs")    x$EFFORT  <- x$EFFORT * 60
                 if(unit == "hours")   x$EFFORT  <- x$EFFORT/60
                 if(unit == "days")    x$EFFORT  <- x$EFFORT / 60 / 24
                 if(unit == "weeks")   x$EFFORT  <- x$EFFORT / 60 / 24 / 7
               }
        return(x)}
        

