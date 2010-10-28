
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# A 'R' ROUTINE FOR THE COUPLING OF VMS AND LOGBOOKS
# WP4 - Lot2 EU tender VMS/LOGBOOKS COUPLING
# author: Francois Bastardie (DTU- Aqua; fba@aqua.dtu.dk)
# January 2010 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#




  
 
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!MERGE LOGBOOKS WITH VMS PER VESSEL!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##


#!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!#
mergeTacsat2EflaloAndDispatchLandingsAtThePingScale <-
           function(logbooks, tacsat, general=list(output.path=file.path("C:"),
                    a.year=2009, visual.check=TRUE, do.wp3=FALSE, speed="segment"), ...){

  lstargs <- list(...)

  #utils--
  an <<- function(x) as.numeric(as.character(x)) # alias
  
  # create required folders for outputs
  cat("if it still doesn't exist, 'results' folder is created in ",general$output.path,"\n")    
  system(paste("mkdir ",file.path(general$output.path),sep=""),intern=TRUE)

  
 #!!!!!!!!!!!!!!!!!!!!!#
 #!!!!!!!!!!!!!!!!!!!!!#
 # utils--
 collapse.all.columns <- function (obj, columns= seq(ncol(obj)) ){
            eval(parse(text=paste('paste(obj[,', paste(columns,collapse='] ,"#", obj[,'), '],sep=\'\')', sep='')))  }
 uncollapse.column <-  function(obj, column="coll"){
            dd<- strsplit(as.character(obj[,column]),"#") ; nco <- length(dd[[1]]) ; dd<- unlist(dd)
            res <- eval(parse(text=paste('data.frame(',paste('dd[seq(',1:nco,',nrow(obj)*nco,by=nco)]', collapse=','),')')))
            colnames(res) <- paste("col",1:nco,sep='')
            return(res)
            }

  #utils--
  # FUNCTION TO CREATE A SPATIAL GRID
  # 'xx' have a 'SI_LATI' and a 'SI_LONG' columns
  assignPointsToSpatialGrid <- function(xx){

    xx <- xx[,!colnames(xx) %in% c("icessquare","icessquare.vms") ]  # remove
    xx <- cbind.data.frame(xx, icessquare= rep(0,nrow(xx)))


    an <- function(x) as.numeric(as.character(x))
    rlong      <- range(an(xx$SI_LONG),na.rm=T)
    vect.long  <- signif(seq(floor(rlong[1]), ceiling(rlong[2]), by=1),4)   # long (x)
    label.long <- rep(paste(rep(LETTERS,each=10),0:9,sep=""),each=1)
    names(label.long) <- signif(seq(-50, 209, by=1),4)   # long (x)
    label.long <- label.long[!is.na(names(label.long))]  # => correspondance long (-50 to 209) / sq letter (A0 to Z9)
    label.long <- label.long[as.character(vect.long)]
    rlat      <- range(an(xx$SI_LATI), na.rm=T)
    vect.lat   <- signif(seq(floor(rlat[1]), ceiling(rlat[2]),by=0.5),4) # lat  (y)
    label.lat  <- rep(paste(seq(1,75,1)),each=1)
    names(label.lat) <-   paste(signif(seq(36,73, by=0.5),4))
    label.lat <- label.lat[!is.na(names(label.lat))] # => correspondance lat (36 to 73) / sq number (1 to 75)
    label.lat <- label.lat[as.character(vect.lat)]
    vect.label <- paste(rep(label.lat,each=length(label.long)),"",label.long,sep="")
    xx[,"SI_RECT"] <- paste(label.lat [findInterval(an(xx[,"SI_LATI"]), vect.lat)] , label.long [findInterval(an(xx[,"SI_LONG"]), vect.long)], sep="")

   return(xx)
   }


   #!!!!!!!!!!!!!!!!!!!!!#
   #utils--
   # for managing NA on logbook side
   # (from vms trip.sq without corresponding logbook trip.sq e.g. because no declaration in sq because only steaming time inside)
   # we need to inform back the specificity of the vessel from logbook using info from the same trip i.e. vesselid+bk.tripnum
   retrieveOnBkSide <- function(merged, type.data){
      idx <- which(merged$LE_MET_level6=="NA")
      merged.NA <- merged[idx,] # input (only the trip.sq with NA for the logbook part)

      for (td in type.data){
         map <- tapply(merged[, td ], paste(merged$VE_REF, merged$bk.tripnum),
                             function(i) {ss<- unique(as.character(i)) ; ss[ss!="NA"][1]})
         merged.NA[, td ] <- factor(paste(merged.NA$VE_REF,merged.NA$bk.tripnum))
         levels(merged.NA[, td ]) <- map[levels(merged.NA[, td ])]
         }
      if(nrow(merged.NA)>0) merged.NA$flag <- 4 # flag on meth
      merged[idx,] <- merged.NA # output
      return(merged)
      }

      #!#!##!#!##!#!##!#!##!#!##!#!#
      #!#!##!#!##!#!##!#!##!#!##!#!#
      #!#!##!#!##!#!##!#!##!#!##!#!#
      #!#!##!#!##!#!##!#!##!#!##!#!#
      #!#!##!#!##!#!##!#!##!#!##!#!#
      all.vesselid     <- as.character(unique(logbooks[an(logbooks$VE_LEN)>=0,]$VE_REF)) 
      all.vesselid     <- all.vesselid[!is.na(all.vesselid)] # e.g. when VE_LEN at NA exists     
      if(length(lstargs$a.vesselid)!=0) all.vesselid <- lstargs$a.vesselid 
       # => IF ARG INFORMED, THEN KEEP ONLY ONE OR SEVERAL VESSELS AS NEEDED....

      for(a.vesselid in all.vesselid){  # PER VESSEL
                cat(paste(a.vesselid,"\n", sep="" ))
       
         #----------
         #----------
         #----------
         #----------
         #----------
         #----------
         # LOGBOOK INPUT
         logbk.this.vessel            <- logbooks[logbooks$VE_REF %in% a.vesselid,]
         logbk.this.vessel$LE_RECT    <- factor(logbk.this.vessel$LE_RECT)
         logbk.this.vessel$VE_REF     <- factor(logbk.this.vessel$VE_REF)

         #   add mid-time and bk.tripnum in eflalo
           # departure time
           ctime <- strptime(  paste(logbk.this.vessel$FT_DDAT, logbk.this.vessel$FT_DTIME) , tz='GTM',  "%e/%m/%Y %H:%M" )
           logbk.this.vessel <- cbind.data.frame(logbk.this.vessel, date.in.R.dep=ctime)
           # arrival time
           ctime <- strptime(  paste(logbk.this.vessel$FT_LDAT, logbk.this.vessel$FT_LTIME) , tz='GTM',  "%e/%m/%Y %H:%M" )
           logbk.this.vessel <- cbind.data.frame(logbk.this.vessel, date.in.R.arr=ctime)
           # catch.date
           ctime <- strptime(  paste(logbk.this.vessel$LE_CDAT) , tz='GTM',  "%e/%m/%Y" )
           logbk.this.vessel <- cbind.data.frame(logbk.this.vessel, date.in.R.cat=ctime)

           # mid time bk trips
           mid.time <- rep(NA, nrow(logbk.this.vessel))
           dep <- logbk.this.vessel$date.in.R.dep +10  # we artificially add +10min because bug in R if mid-time is 00:00:00
           arr <- logbk.this.vessel$date.in.R.arr +1

           for(r in 1:length(dep)){
              mid.time[r] <- as.character(seq(from=dep[r], to=arr[r], length.out = 3)[2])
              }
           logbk.this.vessel$mid.time           <-  mid.time
           logbk.this.vessel$bk.tripnum         <-  factor(mid.time) # init
           levels(logbk.this.vessel$bk.tripnum) <- 1:length(logbk.this.vessel$bk.tripnum) # assign a bk.tripnum code from mid.time
     
           
         #=> LOGBOOK (EFLALO) INPUT REQUIRES AT LEAST,
         #     'VE_REF',  FT_DDAT, FT_DTIME, FT_LDAT, FT_LTIME, FT_CDAT,
         #  'LE_SP_KG' (etc.), 'LE_RECT', 'VE_FLT' AND 'LE_MET_level6', 'LE_GEAR' COLUMNS
         #

         #----------
         #----------
         #----------
         #----------
         #----------
         #----------
         # VMS INPUT: load traj with 'at sea' pings SI_STATE informed
         # ABSOLUTELY REQUIRED: c("VE_REF","SI_LATI","SI_LONG", "SI_DATE", "SI_TIME", "SI_FT", "SI_HARB", "SI_STATE")
    

         if(a.vesselid %in% unique(tacsat$VE_REF)){
      
         tacsat.this.vessel <- tacsat[tacsat$VE_REF == a.vesselid,] # subset for this vessel
         tacsat.this.vessel$VE_REF <- factor(tacsat.this.vessel$VE_REF)
         
 
         # if does not exist, add date.in.R for handling the time in R
         if(!("date.in.R" %in% colnames(tacsat))){
           ctime <- strptime(  paste(tacsat.this.vessel$SI_DATE, tacsat.this.vessel$SI_TIME) , 
                                 tz='GTM',   "%e/%m/%Y %H:%M" )
           tacsat.this.vessel <- cbind.data.frame(tacsat.this.vessel, date.in.R=ctime)
         }

         # keep only the essential
         vms.this.vessel  <- tacsat.this.vessel [, c("VE_REF","SI_LATI","SI_LONG", 
                          "date.in.R","SI_FT", "SI_HARB", "SI_STATE")]
         rm(tacsat.this.vessel); gc(reset=TRUE)                  
         vms.this.vessel$VE_REF   <- factor(vms.this.vessel$VE_REF)

         vms.this.vessel$idx  <- 1:nrow(vms.this.vessel) # label for each ping


        
 
         # filter if vessel with a bad vms
         to.remove.because.deficient.vms <- any(is.na(vms.this.vessel$SI_FT))
         to.remove.because.not.enough.vms.trips <- length(unique(vms.this.vessel$SI_FT))< 2  # nb vms trips < 2
         to.remove.because.pble.lgbk <- length(unique(logbk.this.vessel$FT_REF))< 2  # nb logbk trips < 2
         if(length(unique(vms.this.vessel$SI_FT))<2) warning('need more than 1 trip in SI_FT')
         a.flag <- to.remove.because.deficient.vms ||  to.remove.because.not.enough.vms.trips || to.remove.because.pble.lgbk
         
         ## remove bk.tripnum and mid.time if it exists
         vms.this.vessel <- vms.this.vessel[, !colnames(vms.this.vessel) %in% c("bk.tripnum", "mid.time")]



        if(a.flag==FALSE) {  # i.e. vms-equipped

           if(all(is.na(vms.this.vessel$SI_STATE)) && general$do.wp3==FALSE)
                  stop('the SI_STATE column has to be informed before making the merging')
     
         # alias
         .logbk <- logbk.this.vessel
         .vms   <- vms.this.vessel





         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         #!  DO THE VIRTUAL MERGING  #!#!#!#!#!#!#!#!#!#!#
         #!#!#!#!#!#!#!#!!#!#!#!#!#!#!#!!#!#!#!#!#!#!#!#!#
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#

         #!!!!!!!!!!!!!!!!!!#
         #!!!!!!!!!!!!!!!!!!#
         # -If IT DOES NOT EXIST YET-,
         # FIND THE MID-TIME OF VMS TRIPS
         if(any(colnames(.vms)%in%"date.in.R")){
           if(!any(colnames(.vms)%in%"date.in.R.dep")){
            # find and add the first point of each trip
           .vms$start.trip <- c(1,diff(.vms[,"SI_FT"]))
           .vms$end.trip <- c(diff(.vms[,"SI_FT"]),0)
           .vms[.vms$start.trip>0, "start.trip"] <- .vms[.vms$start.trip>0, "SI_FT"]
           .vms[.vms$end.trip>0, "end.trip"] <- .vms[.vms$end.trip>0, "SI_FT"]

           tmp <- .vms[.vms$start.trip>0,]
           tmp <- tmp[,c("VE_REF","date.in.R","SI_FT")]
           tmp2 <- .vms[.vms$end.trip>0,]
           tmp2 <- tmp2[,c("VE_REF","date.in.R","SI_FT")]
           .vms <- .vms[,!colnames(.vms) %in% c("start.trip", "end.trip")] # remove tool columns
           table.midtime <- merge(tmp, tmp2, by.x="SI_FT", by.y="SI_FT")
           table.midtime <- table.midtime[, c("SI_FT","VE_REF.x","date.in.R.x","date.in.R.y") ]
           colnames(table.midtime) <- c("SI_FT","VE_REF","date.in.R.dep","date.in.R.arr")
           } else{
           table.midtime <- .vms[, c("SI_FT","VE_REF","date.in.R.dep","date.in.R.arr") ]
           table.midtime <- table.midtime[!duplicated(data.frame(table.midtime$SI_FT, table.midtime$VE_REF)),]
           }
         } else{stop("no 'date.in.R' found in vms")}
         mid.time <- rep(0, nrow(table.midtime))
         for(r in 1: nrow(table.midtime)){
           mid.time[r] <-  as.character(seq(from=table.midtime$date.in.R.dep[r], to=table.midtime$date.in.R.arr[r], length.out = 3)[2])
         
         }
         table.midtime$mid.time <-  mid.time
         if(!any(colnames(.vms)%in%"mid.time")){ # here we are...
              .vms <- merge(.vms, table.midtime[,c("SI_FT","mid.time")], by.x="SI_FT", by.y="SI_FT")
         }



        #!!!!!!!!!!!!!!!!!!#
        #!!!!!!!!!!!!!!!!!!#
        # ASSIGN A 'BK.TRIPNUM' FROM LOGBOOK TO EACH VMS TRIP
         trunk <-1 # trunk give the part of the year to be plotted (1 to 5)
         # visual check
         if(general$visual.check){
            windows(width=8, height=4)
            ltrunk <- (nrow(table.midtime)/5)
            idxtrunk <-  (trunk+(trunk-1)*ltrunk):(trunk*ltrunk)
        #    plot(table.midtime$date.in.R.dep[idxtrunk],rep(1,length(table.midtime$date.in.R.dep[idxtrunk])),
             plot(table.midtime$date.in.R.dep,rep(1,length(table.midtime$date.in.R.dep)),
                 ylim=c(0,0.52), type="n", ylab="", axes=FALSE)
            r <- as.POSIXct(round(range(table.midtime$date.in.R.dep), "days"))
            axis.POSIXct(1, at=seq(r[1], r[2], by="month"), format="%e%b%y:%H:%M")
            axis(2, at=c(0.5,0.1),labels=c("VMS","LOGBOOK"))

            for(i in 1:nrow(table.midtime))  {
              segments(as.POSIXct(table.midtime$date.in.R.dep[i]), 0.5, as.POSIXct(table.midtime$date.in.R.arr[i]), 0.5, col=1)
              points(as.POSIXct(table.midtime$mid.time[i]), 0.5, col=1)
              text(as.POSIXct(table.midtime$mid.time[i]), 0.52, table.midtime$SI_FT[i], cex=0.5, col=1)
    
            }
            tmp <- .logbk[, c("date.in.R.dep","date.in.R.arr", "mid.time", "bk.tripnum")]
            tmp <- tmp[!duplicated(tmp$mid.time), ]
            for(i in 1:nrow(tmp)){
              segments(as.POSIXct(tmp$date.in.R.dep[i]), 0.1, as.POSIXct(tmp$date.in.R.arr[i]), 0.1, col=1)
              points(as.POSIXct(tmp$mid.time[i]), 0.1, col=1)
              text(as.POSIXct(tmp$mid.time[i]), 0.0785, tmp$bk.tripnum[i], cex=0.5, col=1)
            }
          }
         

          # THE CORE CODE: compare bk$mid.time and vms$mid.time
          # find the nearest bk$mid.time for each vms$mid.time
          # and then change levels
          # (so, for each mid.time in vms, a bk.tripnum will be find)
          # (so, no lines in vms without a bk.tripnum from bk...)
          fa1 <- levels(factor(.vms$mid.time))
          new.levels <- fa1
          fa2 <-  levels(factor(.logbk$mid.time))
          for(i in 1:length(fa1)) { # for each level in vms
             tmp <-  abs(as.numeric( as.POSIXct(fa2) - as.POSIXct(fa1)[i] ))
             if(all(is.na(tmp))) tmp <- abs(as.numeric( as.Date(fa2) - as.Date(fa1)[i] )) # debug the R bug in case of mid-time at 00:00 hour
             new.levels[i] <- fa2 [which(tmp == min(tmp, na.rm=T) )]  # find the nearest level in logbook
          }
          .vms$mid.time <- factor(as.character(.vms$mid.time))
          sauv <- .vms$mid.time
          levels(.vms$mid.time) <- new.levels # and change mid.time in vms to force the merging

          # finally, replace levels by the bk.tripnum
          tmp <-  .logbk[.logbk$mid.time %in% .vms$mid.time , c("bk.tripnum","mid.time")]
          tmp2 <- tmp[!duplicated(tmp$bk.tripnum),]
          idx <- match(levels(.vms$mid.time),tmp2$mid.time)
          .vms$bk.tripnum <- .vms$mid.time # init
          levels(.vms$bk.tripnum) <- as.character(tmp2$bk.tripnum )   [idx]


          if(general$visual.check){
            for(i in 1: nrow(.vms))  {
               arrows(as.POSIXct( sauv[i]), 0.5 ,as.POSIXct( .vms$mid.time[i]),0.1, length=0.1)
            }
          }

          if(general$visual.check){
            ve <- as.character(.logbk$VE_REF[1])
            savePlot(filename = file.path(general$output.path,
                            paste("assign_eflalo_tripnum_to_vms_",ve,"_",general$a.year,".jpeg",sep="")),type ="jpeg")
           dev.off()
          }

     
        ## ADD A WARNING IN CASE OF LONG (UNREALISTIC) TRIPS ##
        diff.date <- table.midtime$date.in.R.arr - table.midtime$date.in.R.dep    # if at least one trip >30 days
        if(attributes(diff.date)$units=="secs")  idx <- which((((diff.date)/3600)/24) >30)  
        if(attributes(diff.date)$units=="hours")  idx <- which((((diff.date)/1)/24) >30)  
        attributes((table.midtime$date.in.R.arr - table.midtime$date.in.R.dep ))
        if (length( idx) >0){
             cat(paste("at least one vms trip > 30 days detected! check harbours...", "\n", sep=""))
            suspicious <- .vms[.vms$SI_FT %in%  table.midtime$SI_FT[idx] ,]
            tmp <- table(suspicious$SI_LATI)
            lat.suspicious <- names(tmp[tmp>5]) 
            if(length(lat.suspicious)!=0) cat(paste("potential harbour likely near lat ",lat.suspicious,"\n",sep=""))
            tmp <- table(suspicious$SI_LONG)
            long.suspicious <- names(tmp[tmp>5]) 
            if(length(long.suspicious)!=0) cat(paste("potential harbour likely near long ",long.suspicious,"\n",sep=""))
            }  # if at least one trip >30 days
        rm(table.midtime) ; gc(reset=TRUE)  

     
       
         .logbk$mid.time    <- factor(.logbk$mid.time)
         .logbk$bk.tripnum  <- factor(.logbk$bk.tripnum)
         .vms$mid.time      <- factor(.vms$mid.time)
         .vms$bk.tripnum    <- factor(.vms$bk.tripnum)

         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         #! ASSIGN A 'SI_FT' FROM VMS TRIP NUM TO  #!#!#!#
         #! LOGBOOK TRIPS WITH NO VMS CORRESPONDANCE #!#!#
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
       dep.bk.not.in.vms  <- unique( .logbk$bk.tripnum [ !( .logbk$bk.tripnum %in%  .vms$bk.tripnum )  ] )
       if(length(dep.bk.not.in.vms)!=0){
         # bk tripnum from dep not in vms
         idx <-  .logbk$bk.tripnum %in% dep.bk.not.in.vms
         bk  <- .logbk[idx,] [order(.logbk[idx,]$date.in.R.dep),]
         if(!"mid.time" %in% colnames(.vms)){
            vms <- .vms  [order(.vms$date.in.R.dep),]
            mid.time <- rep(NA, nrow(vms))
            for(r in 1: nrow(vms)){
              mid.time[r] <-  as.character(seq(from=vms$date.in.R.dep[r], to=vms$date.in.R.arr[r], length.out = 3)[2])
            }
            vms$mid.time <-  mid.time
         } else{ vms <- .vms[order(.vms$mid.time),]}
         #1- compare bk$mid.time and vms$mid.time
         # find the nearest vms$mid.time for each bk$mid.time
         # and then change levels
         # (so for each mid.time in bk, a tripnum will be find)
         # (so no lines in bk without a tripnum...)
         fa1 <- levels(factor(bk$mid.time))
         new.levels <- fa1
         fa2 <-  levels(factor(vms$mid.time))
         for(i in 1:length(fa1)) { # for each level in logbk
          tmp <-  abs(as.numeric( as.POSIXct(fa2) - as.POSIXct(fa1)[i] ))
          new.levels[i] <- fa2 [which(tmp == min(tmp, na.rm=T) )]  # find the nearest level in vms
         }
         bk$mid.time <- factor(as.character(bk$mid.time))
         levels(bk$mid.time) <- new.levels # and change mid.time in logbk to force the merging

         # finally, replace levels by the tripnum
         # (note: a same bk.tripnum in vms can have different mid.time
         # due to the first merging of vms to logbk in the vms analysis)
         tmp <-  vms[vms$mid.time %in% bk$mid.time , c("bk.tripnum","mid.time")]
         tmp2 <- tmp[!duplicated(data.frame(tmp$bk.tripnum,tmp$mid.time)),]
         idx2 <- match(levels(bk$mid.time), tmp2$mid.time)
         bk$bk.tripnum <- bk$mid.time # init
         levels(bk$bk.tripnum)  <- as.character(tmp2$bk.tripnum) [idx2]

         # output
         bk$mid.time   <- as.character(bk$mid.time)
         bk$bk.tripnum <- as.character(bk$bk.tripnum)
         .logbk[idx,][order(.logbk[idx,]$date.in.R.dep),]  <- bk
         }




         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         # ASSIGN A RECTANGLE TO EACH PING #!#!#!#!#!#!#!#
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         .vms   <- assignPointsToSpatialGrid(xx=.vms)
        
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         # COMPUTE EFFORT.MINS      !#!#!#!#!#!#!#!#!#!#!#
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
          .vms <- .vms[order(.vms$date.in.R),]
          .vms$LE_EFF_VMS <- abs(c(0, as.numeric(.vms[-nrow(.vms),"date.in.R"] - 
                                        .vms[-1,"date.in.R"], units="mins")))
           start.trip <- c(1,diff(.vms[,"SI_FT"]))
          .vms[start.trip!=0, "LE_EFF_VMS"] <- 0  # just correct for the trip change points


         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         # ASSIGN FISHING/NON-FISHING (optional)!#!#!#!#!#
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         if(general$do.wp3 && general$speed=="segment")
            .vms <- segmentSpeedTacsat (tacsat=.vms, vessels=a.vesselid, 
                                  force.lower.bound=0.5, general=list(
                                   output.path=general$output.path,visual.check=TRUE))
                #=> (semi)automatic detection of the fishing peak
                # (put here because the LE_GEAR need to be informed)
         # alternative TO DO:
         #if(general$do.wp3 && general$speed=="lookuptable")
         #   .vms <- lookupSpeedTacsat (tacsat=.vms, vessels=a.vesselid)



         rm(er); rm(xx) ; gc(reset=TRUE)
                                           
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         # SET UP PRIMARY KEYS FOR MERGING!#!#!#!#!#!#!#!#
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         .logbk$bk.tripnum <- factor(.logbk$bk.tripnum )
         .logbk$bk.tripnum.sq <- paste(.logbk$bk.tripnum, ".", .logbk$LE_RECT, sep='') # caution:redefine
         .logbk$bk.tripnum.sq.day <- paste(.logbk$bk.tripnum, ".", .logbk$LE_RECT,".",.logbk$date.in.R.cat, sep='') # caution:redefine
         .vms$bk.tripnum <- factor(.vms$bk.tripnum)
         .vms$bk.tripnum.sq <- paste(.vms$bk.tripnum, ".", .vms$SI_RECT, sep='') # caution:redefine
         .vms$bk.tripnum.sq.day <- paste(.vms$bk.tripnum, ".", .vms$SI_RECT,".", format(.vms$date.in.R,  '%Y-%m-%d'), sep='') # caution:redefine

         # for gear, if several gears inside a same trip,
         #  it is problematic because we have to assume a split of total effort or toal nb of ping between gears...


           #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
           # AGGREGATE WEIGHT PER SPECIES !#!#!#!#!#!#!#!#!#
           #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
           nm         <- names(.logbk)
           idx.col.w  <- grep('KG', nm) # index columns with species weight
           idx.col.v  <- grep('EURO', nm) # index columns with species value
           idx.col    <- c(idx.col.w, idx.col.v)
             # AGGREGATE WEIGHT (OR VALUE) PER SPECIES PER BK.TRIPNUM
              agg.logbk.this.vessel.method.1  <- aggregate(.logbk[,idx.col],
                      list(.logbk$bk.tripnum, 
                              .logbk$VE_REF, .logbk$VE_KW, .logbk$VE_FLT,  .logbk$LE_MET_level6, .logbk$LE_GEAR), sum, na.rm=TRUE )
              colnames(agg.logbk.this.vessel.method.1) <- 
                           c("bk.tripnum", "VE_REF", "VE_KW", "VE_FLT", "LE_MET_level6","LE_GEAR", nm[idx.col] )
             # AGGREGATE WEIGHT (OR VALUE) PER SPECIES PER BK.TRIPNUM.SQ
              agg.logbk.this.vessel.method.2  <- aggregate(.logbk[,idx.col],
                      list(.logbk$bk.tripnum.sq, 
                              .logbk$VE_REF, .logbk$VE_KW, .logbk$VE_FLT,  .logbk$LE_MET_level6, .logbk$LE_GEAR), sum, na.rm=TRUE )
              colnames(agg.logbk.this.vessel.method.2) <- 
                           c("bk.tripnum.sq", "VE_REF", "VE_KW", "VE_FLT","LE_MET_level6" ,"LE_GEAR", nm[idx.col])
             # AGGREGATE WEIGHT (OR VALUE) PER SPECIES PER BK.TRIPNUM.SQ.DAY (NOTE: SO, 'LE_SEQNUM' IS AGGREGATED HERE)
              agg.logbk.this.vessel.method.3  <- aggregate(.logbk[,idx.col],
                      list(.logbk$bk.tripnum.sq.day, 
                             .logbk$VE_REF, .logbk$VE_KW, .logbk$VE_FLT,  .logbk$LE_MET_level6, .logbk$LE_GEAR), sum, na.rm=TRUE )
              colnames(agg.logbk.this.vessel.method.3) <- 
                          c("bk.tripnum.sq.day", "VE_REF", "VE_KW", "VE_FLT","LE_MET_level6","LE_GEAR",  nm[idx.col])


             #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
             # MERGING WITH VMS PER TRIP !!!!!!!!!!#!#!#!#!#!#
             #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
             do.merging <- function(method="bk.tripnum", .logbk, .vms, general){


              # IF BY PING-------------
              # find total nb of FISHING ping per tripnum from vms  # used for method 1  'bk.tripnum'
              if(method=="bk.tripnum"){
               .vms$count.fping.trip  <- factor(.vms$bk.tripnum)  # init
              count.fping.trip <- table(.vms[.vms$SI_STATE==1,]$bk.tripnum)
              # => COUNT nb of FISHING pings per bk.tripnum because each weight will be repeated by ping after merging
              levels(.vms$count.fping.trip) <- count.fping.trip[levels(.vms$count.fping.trip)]  # mapping
              .vms[.vms$SI_STATE==2,]$count.fping.trip <- NA
              # => COUNT nb of gears per bk.tripnum because each ping will be repeated by gear after merging
              count.gr.trip <- tapply(.logbk$LE_GEAR, .logbk$bk.tripnum, function(x) length(unique(x)))
              .logbk$count.gr.trip <- count.gr.trip[.logbk$bk.tripnum]  # mapping
              }


              # find total nb of FISHING ping per trip-icessquare from vms  # used for method 2   'bk.tripnum.sq'
              if(method=="bk.tripnum.sq"){
              .vms$count.fping.trip.sq  <- factor(.vms$bk.tripnum.sq)  # init
              count.fping.trip.sq <- table(.vms[.vms$SI_STATE==1,]$bk.tripnum.sq) # COUNT nb of FISHING pings per bk.tripnum.sq
              levels(.vms$count.fping.trip.sq) <- count.fping.trip.sq[levels(.vms$count.fping.trip.sq)]  # mapping
              if(any('2' %in% unique(.vms$SI_STATE))) .vms[.vms$SI_STATE==2,]$count.fping.trip.sq <- NA
              # => COUNT nb of gears per bk.tripnum.sq because each ping will be repeated by gear after merging
              count.gr.trip.sq <- tapply(.logbk$LE_GEAR, .logbk$bk.tripnum.sq, function(x) length(unique(x)))
              .logbk$count.gr.trip.sq <- count.gr.trip.sq[.logbk$bk.tripnum.sq]  # mapping
              }

              # find total nb of FISHING ping per trip-icessquare-day from vms  # used for method 3   'bk.tripnum.sq.day'
              if(method=="bk.tripnum.sq.day"){
              .vms$count.fping.trip.sq.day  <- factor(.vms$bk.tripnum.sq.day)  # init
              count.fping.trip.sq.day <- table(.vms[.vms$SI_STATE==1,]$bk.tripnum.sq.day) # COUNT nb of FISHING pings per bk.tripnum.sq.day
              levels(.vms$count.fping.trip.sq.day) <- count.fping.trip.sq.day[levels(.vms$count.fping.trip.sq.day)]  # mapping
              if(any('2' %in% unique(.vms$SI_STATE))) .vms[.vms$SI_STATE==2,]$count.fping.trip.sq.day <- NA
              # => COUNT nb of gears per bk.tripnum.sq.day because each ping will be repeated by gear after merging
              count.gr.trip.sq.day <- tapply(.logbk$LE_GEAR, .logbk$bk.tripnum.sq.day, function(x) length(unique(x)))
              .logbk$count.gr.trip.sq.day <- count.gr.trip.sq.day[.logbk$bk.tripnum.sq.day]  # mapping}
              }



              # do the merging between .logbk and .vms according to
              #  meth1: 'bk.tripnum' OR meth2: 'bk.tripnum.sq' OR meth3: 'bk.tripnum.sq.day'
              # need to use a trick to avoid "out of memory" doing the merge()
              coln.idx1 <- which(!colnames(.logbk)%in%c("VE_REF", method))
              coln1 <- colnames(.logbk)[coln.idx1]
              tmp1 <- data.frame(coll= collapse.all.columns  (.logbk, columns= coln.idx1  ),
                         VE_REF=.logbk$VE_REF, a.method= .logbk[,method] ) #.logbk
              coln.idx2 <- which(!colnames(.vms)%in%c("VE_REF", method))
              coln2 <- colnames(.vms)[coln.idx2]
              tmp2 <- data.frame(coll2= collapse.all.columns  (.vms, columns=  coln.idx2 ),
                         VE_REF=.vms$VE_REF, a.method= .vms[,method] )  #.vms
              tmp1[,"a.method"] <- factor(tmp1[,"a.method"] )
              tmp2[,"a.method"] <- factor(tmp2[,"a.method"] )

              merged.this.vessel <- merge(tmp1, tmp2, all.x=TRUE, all.y=TRUE, suffixes = c(".bk",".vms"))
              #=> so, with all.y = TRUE, the vms records without corresponding logbk records are kept and NA are produced on the logbook part
              #=> so, with all.x = TRUE, the logbk records  without corresponding vms records are kept and NA are produced on the vms part
              merged.this.vessel$coll <- replace(as.character(merged.this.vessel$coll),is.na(merged.this.vessel$coll), paste(rep("NA",length(coln1)),collapse="#"))
              merged.this.vessel$coll <- factor(merged.this.vessel$coll)
              #=> adapt 'coll' to get a vector of NA (NA in case of 'in vms but not in logbook')
              merged.this.vessel$coll2 <- replace(as.character(merged.this.vessel$coll2),is.na(merged.this.vessel$coll2), paste(rep("NA",length(coln2)),collapse="#"))
              # adapt 'coll2' to get a vector of NA (NA in case of 'in logbook but not in vms')
              merged.this.vessel$coll2 <- factor(merged.this.vessel$coll2)
              colnames(merged.this.vessel)[colnames(merged.this.vessel)%in%"a.method"] <- method

              tmp3 <- uncollapse.column(merged.this.vessel, column="coll")  # logbk
              tmp4 <- uncollapse.column(merged.this.vessel, column="coll2") # vms
              tmp5 <- cbind.data.frame(merged.this.vessel[,c("VE_REF", method)], tmp3, tmp4)
              colnames(tmp5) <- c("VE_REF", method, coln1, coln2)
              merged.this.vessel <- tmp5

              # note: at this stage some few lgbk records could have got 
              # no correpondance in vms (misreporting of area)=> NA on vms side
              # and then the landing weight of these records are possibly lost because count.ping at NA...
              # so we can choose to correct (see ** below) to keep these land. weight
              # the remaining loss in weight will come from the matching records having catches but
              # without fishing pings (i.e. only steaming pings)!

          

              if(FALSE){
              # conservation of catches?
              # detect possible weight landed while no feffort detected from vms
                   # find bk.tripnum with some NA
                   vv<- an(unique(merged.this.vessel[merged.this.vessel$count.fping.trip=="NA","bk.tripnum"]))
                   # then, find bk.tripnum with at least one no NA
                   no.vv<- an(unique(merged.this.vessel[merged.this.vessel$count.fping.trip!="NA","bk.tripnum"]))
                   tripnum.all.na.inside <- vv[!vv%in%no.vv] # trip num without at least one count.fping!
                   # so, deduce loss in weight
                   zz<- merged.this.vessel[merged.this.vessel$bk.tripnum %in% tripnum.all.na.inside,]
                   sum(an(unique(zz$LE_KG_COD)), na.rm=TRUE)

                cat(paste("weight loss for ", general$sp.to.keep[1]," (vms failure in fishing/steaming detection): ",
                      sum(an(unique(zz$LE_KG_COD)), na.rm=TRUE),"\n", sep="" ))
              
             }  # TO DO**: assign landings to the mid point of the trip for trips with all na inside (i.e. only steaming detected while declared landings) 
                 # (i.e. assign 1 to in count.fping.trip for the mid point)

              # apply the catches re-distribution
              # method 1, 2 and 3: per ping
              # PER PING:
              # ASSUMING EQUAL ALLOCATION BETWEEN FISHING PINGS AND GEARS USE INSIDE A SAME TRIP
              an        <- function(x) as.numeric(as.character(x))
              nm        <- names(merged.this.vessel)
              idx.col.w <- grep('KG', nm) # index columns with species weight
              idx.col.v <- grep('EURO', nm) # index columns with species value
              idx.col <- c(idx.col.w, idx.col.v)
              if(method=="bk.tripnum.sq.day"){
                  merged.this.vessel[merged.this.vessel$SI_LATI=='NA', "count.fping.trip.sq.day"] <- 1 # **correct for loss if in lgk but not in vms
                             merged.this.vessel[,idx.col] <- (apply(merged.this.vessel[,idx.col],2,an) /
                                                        an(merged.this.vessel$count.fping.trip.sq.day)) /
                                                                        an(merged.this.vessel$count.gr.trip.sq.day)
              }
              if(method=="bk.tripnum.sq"){
                  merged.this.vessel[merged.this.vessel$SI_LATI=='NA', "count.fping.trip.sq"] <- 1 # **correct for loss if in lgk but not in vms
                             merged.this.vessel[,idx.col] <- (apply(merged.this.vessel[,idx.col],2,an) /
                                                        an(merged.this.vessel$count.fping.trip.sq)) /
                                                                        an(merged.this.vessel$count.gr.trip.sq)
              }
              if(method=="bk.tripnum"){
                 merged.this.vessel[merged.this.vessel$SI_LATI=='NA', "count.fping.trip"] <- 1 # **correct for loss if in lgk but not in vms
                         # maybe do more by adding unallocated landings to the midpoint of the trip**
                             merged.this.vessel[,idx.col] <- (apply(merged.this.vessel[,idx.col],2,an) /
                                                           an(merged.this.vessel$count.fping.trip) ) /
                                                                        an(merged.this.vessel$count.gr.trip)
              }

      

              return(merged.this.vessel)
              }



             #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
             # MERGING PROCEDURE CHOICE !#!#!#!#!#!#!#!#!#!#!#
             #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
          
                 .logbk   <- agg.logbk.this.vessel.method.3
                 my.split <- function(obj,a.sep="\\.",idx=1) unlist(lapply(strsplit(obj, a.sep),function(x)x[idx]))
                 # reduce the level
                 .logbk$bk.tripnum.sq  <-  paste(my.split(as.character(.logbk$bk.tripnum.sq.day),a.sep="\\.",idx=1),
                                                 my.split(as.character(.logbk$bk.tripnum.sq.day),a.sep="\\.",idx=2),sep='.')
                 # reduce the level
                 .logbk$bk.tripnum     <-        my.split(as.character(.logbk$bk.tripnum.sq),a.sep="\\.",idx=1)
                 # verbose & export
                 # find common keys
                 tripnum.sq.day.logbk            <- .logbk$bk.tripnum.sq.day
                 tripnum.sq.day.vms              <- .vms$bk.tripnum.sq.day
                 tripnum.sq.logbk                <- .logbk$bk.tripnum.sq
                 tripnum.sq.vms                  <- .vms$bk.tripnum.sq
                 tripnum.sq.day.in.vms.and.in.bk <- tripnum.sq.day.vms [tripnum.sq.day.vms %in% tripnum.sq.day.logbk]
                 tripnum.sq.in.vms.and.in.bk     <- tripnum.sq.vms [tripnum.sq.vms %in% tripnum.sq.logbk]
                 .vms.in.bk                      <- .vms[ .vms$bk.tripnum.sq.day %in%  tripnum.sq.day.in.vms.and.in.bk,]
                 .vms.in.bk2                     <- .vms[ .vms$bk.tripnum.sq %in%  tripnum.sq.in.vms.and.in.bk,]
                 in.bk.and.feffort.not.at.0   <- unique(.vms.in.bk[.vms.in.bk$SI_STATE==1,]$bk.tripnum.sq.day)
                 in.bk2.and.feffort.not.at.0   <- unique(.vms.in.bk2[.vms.in.bk2$SI_STATE==1,]$bk.tripnum.sq)
                 
                     # split .vms and .logbk in three blocks
                  # vms with good match => go to meth3
                 .vms.for.meth3         <- .vms [.vms$bk.tripnum.sq.day %in%   in.bk.and.feffort.not.at.0, ]
                  # vms with intermediate match => go to meth2
                 .vms.for.meth2         <- .vms [!(.vms$bk.tripnum.sq.day  %in%   in.bk.and.feffort.not.at.0) &
                                                      (.vms$bk.tripnum.sq    %in%   in.bk2.and.feffort.not.at.0), ]
                  # vms with bad match => go to meth1
                 .vms.for.meth1         <- .vms [!(.vms$bk.tripnum.sq.day  %in%   in.bk2.and.feffort.not.at.0) &
                                                      !(.vms$bk.tripnum.sq  %in%   in.bk2.and.feffort.not.at.0), ]
                  # logbk with good match => go to meth3
                 .logbk.for.meth3       <- .logbk [.logbk$bk.tripnum.sq.day %in%  in.bk.and.feffort.not.at.0, ]
                  # logbk with intermediate match => go to meth2
                 .logbk.for.meth2       <- .logbk [!(.logbk$bk.tripnum.sq.day %in%   in.bk.and.feffort.not.at.0) &
                                                       (.logbk$bk.tripnum.sq %in%  in.bk2.and.feffort.not.at.0), ]
                  # logbk with bad match => go to meth1
                 .logbk.for.meth1       <- .logbk [!(.logbk$bk.tripnum.sq.day %in%   in.bk.and.feffort.not.at.0) &
                                                       !(.logbk$bk.tripnum.sq %in%  in.bk2.and.feffort.not.at.0), ]

                 suppressWarnings(rm(merged1, merged2, merged3)) # clear
                 #!! METH1 !!#
                 if(nrow(.logbk.for.meth1)!=0 && nrow(.vms.for.meth1)!=0 ) {
                    # remove useless cols and aggregate according to the key 'bk.tripnum'
                    .logbk.for.meth1 <- .logbk.for.meth1[, !colnames(.logbk.for.meth1)%in% c("bk.tripnum.sq.day","bk.tripnum.sq")]
                    nm        <- names(.logbk.for.meth1)
                    idx.col.w <- grep('KG', nm) # index columns with species weight
                    idx.col.v <- grep('EURO', nm) # index columns with species value
                    idx.col <- c(idx.col.w, idx.col.v)
                    .logbk.for.meth1   <- aggregate(.logbk.for.meth1 [,idx.col],
                                 list(.logbk.for.meth1$VE_REF, .logbk.for.meth1$bk.tripnum,
                                            .logbk.for.meth1$VE_FLT, .logbk.for.meth1$VE_KW, .logbk.for.meth1$LE_MET_level6, .logbk.for.meth1$LE_GEAR), sum, na.rm=TRUE)
                    colnames(.logbk.for.meth1) <- c("VE_REF","VE_KW", "bk.tripnum", "VE_FLT", "LE_MET_level6", "LE_GEAR", nm[idx.col])
                    # do.merging
                    merged1  <- do.merging(method="bk.tripnum", .logbk.for.meth1, .vms.for.meth1, general)
                    # add meth flag
                    merged1$flag <- 1 # meth 1
                    }
                 #!! METH2 !!#
                 if(nrow(.logbk.for.meth2)!=0 && nrow(.vms.for.meth2)!=0 ) {
                    # remove useless cols and aggregate according to the key 'bk.tripnum.sq'
                    .logbk.for.meth2 <- .logbk.for.meth2[, !colnames(.logbk.for.meth2)%in% c("bk.tripnum.sq.day","bk.tripnum")]
                    nm        <- names(.logbk.for.meth2)
                    idx.col.w <- grep('KG', nm) # index columns with species weight
                    idx.col.v <- grep('EURO', nm) # index columns with species value
                    idx.col <- c(idx.col.w, idx.col.v)
                    .logbk.for.meth2   <- aggregate(.logbk.for.meth2 [,idx.col],
                                 list(.logbk.for.meth2$VE_REF, .logbk.for.meth2$VE_KW, .logbk.for.meth2$bk.tripnum.sq,
                                            .logbk.for.meth2$VE_FLT, .logbk.for.meth2$LE_MET_level6, .logbk.for.meth2$LE_GEAR), sum, na.rm=TRUE)
                    colnames(.logbk.for.meth2) <- c("VE_REF", "VE_KW", "bk.tripnum.sq", "VE_FLT", "LE_MET_level6", "LE_GEAR", nm[idx.col])
                    # do.merging
                    merged2 <- do.merging(method="bk.tripnum.sq", .logbk.for.meth2, .vms.for.meth2, general)
                    # add meth flag
                    merged2$flag <- 2 # meth 2
                 }
                 #!! METH3 !!#
                 if(nrow(.logbk.for.meth3)!=0 && nrow(.vms.for.meth3)!=0 ) {
                    # do.merging
                    merged3 <- do.merging(method="bk.tripnum.sq.day", .logbk.for.meth3, .vms.for.meth3, general)
                   # add meth flag
                    merged3$flag <- 3 # meth 3
                 }

                 # bind the three blocks
                 merged <- NULL ; colnm <- NULL
                 for(i in 1: 3){
                   a.table <- try(get(paste('merged',i,sep='')), silent=TRUE)
                   if(class(a.table)!="try-error"){
                     a.table <- a.table[, !colnames(a.table) %in%
                                  c("count.fping.trip.sq.day","count.fping.trip.sq","count.fping.trip",
                                      "tot.fish.effort.trip","tot.fish.effort.trip.sq",
                                         "count.gr.trip", "count.gr.trip.sq", "count.gr.trip.sq.day")] # remove tool columns
                     if(i==1) colnm <-  colnames(a.table) ; if(is.null(colnm)) colnm <-  colnames(a.table)
                     merged <- rbind.data.frame (merged, a.table[, colnm])
                     }
                   }
                 # if still 'not merging' part, retrieve on NA side i.e. occurs when pings in vms but not in bk
                   merged <- retrieveOnBkSide(merged, type.data=c( "VE_FLT","VE_KW","LE_MET_level6"))  # i.e. when metier=='NA'

     
        # clean up
        rm(a.table, merged1, merged2, merged3, merged.this.vessel,.vms, .logbk, logbk.this.vessel, vms.this.vessel)
        gc(reset=TRUE)
   
        # restore eflalo names
        names(merged)  [names(merged) %in% "bk.tripnum"] <- "FT_REF"


        # restore tacsat names               "%e/%m/%Y %H:%M"
        idx <- merged$date.in.R!='NA' # NA is possible when bk not in vms because bk.tripnum vms may belong to another block than block1
        merged$date.in.R <- as.character( merged$date.in.R)
        merged$date.in.R.date  <- NA
        merged[idx,]$date.in.R.date <- paste(substr(merged[idx,]$date.in.R ,9,10),"/",
                                      substr(merged[idx,]$date.in.R , 6,7), "/", substr(merged[idx,]$date.in.R ,1,4), sep='')
        merged$date.in.R.time  <- NA
        merged[idx,]$date.in.R.time <- paste(substr(merged[idx,]$date.in.R , 12,13),":",
                                      substr(merged[idx,]$date.in.R , 15,16), sep='')
        names(merged)  [names(merged) %in% "date.in.R.date"] <- "SI_DATE"
        names(merged)  [names(merged) %in% "date.in.R.time"] <- "SI_TIME"
    
        # last calculation 
        merged$KW_HOURS <- an(merged$VE_KW) * an(merged$LE_EFF_VMS)
    
        # last clean up 
        merged <- merged[, !colnames(merged) %in% c('idx', 'icessquare')]
    
       # save------------
       save("merged",   file=file.path(general$output.path,
             paste("merged_",  a.vesselid,"_",general$a.year,".RData", sep='')))
       cat(paste("save 'merged'...OK\n\n",sep=""))
     

               }else{  # end 'a.flag'
     cat(paste("failure for",a.vesselid,"(probably not vms-equipped)\n"))
     # because no vms for this vessel...
     # TO DO: the logbk way
     #...
       }
     }else{  # end try-error
     cat(paste("failure for",a.vesselid,"(probably not vms-equipped)\n"))
     # because no vms for this vessel...
     # TO DO: the logbk way
     #...
     }




     } # end a.vesselid






return()
}






  ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
  ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
  ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
  ##!!!!!MAIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
  ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
  ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
  ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
 if(FALSE) {


  #\dontrun{
  data(eflalo2)
  data(tacsat)
  data(euharbours)
  # add missing harbours? (still to be fix...)
  #euharbours <- c(euharbours, list(a.harbour1=data.frame(lon='10',lat='10')))
  #euharbours <- c(euharbours, list(a.harbour2=data.frame(,lon='1',lat='1')))

 
  library(doBy)
  tacsat$SI_HARB <- NA
  inHarb <- pointInHarbour(lon=tacsat$SI_LONG,lat=tacsat$SI_LATI,harbours=euharbours,30)
  tacsat$SI_FT <- 1 # init
  idx <- which(inHarb==0)
  tacsat[idx,"SI_FT"] <- cumsum(inHarb) [idx] # add a SI_FT index
  tacsat <- tacsat[which(inHarb==0),] # keep out of harbour points only
  tacsat$SI_STATE <- 2 # init (1: fishing; 2: steaming)
  tacsat$SI_STATE [(tacsat$SI_SP>4 & tacsat$SI_SP<8)] <-1 # fake speed rule for fishing state


                       
  
  # reduce the size of the eflalo data by merging species (e.g. <1 millions euros)
  eflalo <- mergeEflaloSpecies (eflalo2, threshold=1e6) 
  
  # debug
  eflalo2 <- eflalo2[!eflalo2$VE_REF=="NA" &!is.na(eflalo2$VE_REF),]
  
  # TEST FOR A GIVEN SET OF VESSELS
  mergeTacsat2EflaloAndDispatchLandingsAtThePingScale (logbooks=eflalo2, tacsat=tacsat, a.vesselid=c("35", "1518"),
                                                             general=list(output.path=file.path("C:","output"),
                                                                            a.year=2009, visual.check=TRUE,
                                                                             do.wp3=FALSE, speed="segment"))
  # ...OR APPLY FOR ALL VESSELS IN eflalo2
  mergeTacsat2EflaloAndDispatchLandingsAtThePingScale (logbooks=eflalo2, tacsat=tacsat,
                                                             general=list(output.path=file.path("C:","output"),
                                                                            a.year=2009, visual.check=TRUE,
                                                                             do.wp3=FALSE, speed="segment"))
  gc(reset=TRUE)

  # load the merged output table for one vessel
  load(file.path("C:","output","merged_35_2009.RData"))
  
  # ...or bind all vessels
  tmp <- bindAllMergedTables (vessels=c("35", "1518"), species.to.merge=character(), what=character(), 
                      folder = file.path("C:","output"))
 
   # ...and load the merged output table for all vessels
  load(file.path("C:","output","all_merged_2009.RData"))
             
  # map landing of sole from all studied vessels
  df1<- all.merged[,colnames(all.merged)%in% c("SI_LATI","SI_LONG","LE_KG_SOL")]
  df1$SI_LONG <-as.numeric(as.character(df1$SI_LONG))
  df1$SI_LATI <-as.numeric(as.character(df1$SI_LATI))
  vmsGridCreate(df1,nameLon="SI_LONG",nameLat="SI_LATI",cellsizeX =0.05,cellsizeY =0.05)

  # remove steaming points before gridding!
  df2<-df1[-which(is.na(df1$LE_KG_SOL)),]
  vmsGridCreate(df2,nameLon="SI_LONG",nameLat="SI_LATI",cellsizeX =0.05,cellsizeY =0.05)


  # CONVERT TO FISHFRAME FORMAT (might take some time running)
  ff <- mergedTable2Fishframe (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="NLD") )

  #}
 
} # end main
