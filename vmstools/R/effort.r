effort <- function(x,level="trip",unit="hours"){

              #-Add if necessary a datim column
              if(!length(grep("SI_DATIM",colnames(x)))>0){
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
                spl   <- split(x,f=af(x$ID))
                eff   <- lapply(spl,function(y){return(an(difftime(last(y$SI_DATIM),y$SI_DATIM[1],units=unit)))})
                res   <- unlist(eff)
                res   <- data.frame(cbind(names(res),res)); colnames(res) <- c("ID","effort")
              }
              if(all(c("VE_FLT","VE_KW") %in% colnames(x))){
                x$SI_DATIM2 <- as.POSIXct(paste(x$FT_LDAT,  x$FT_LTIME,    sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
                res   <- an(difftime(x$SI_DATIM2,x$SI_DATIM,units=unit))
                res   <- data.frame(cbind(x$ID,res)); colnames(res) <- c("ID","effort")
              }

              res2    <- merge(x,res,by="ID",all=T)
        return(res2)}