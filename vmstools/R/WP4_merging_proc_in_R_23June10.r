
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


#!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!#
#utils--
# FUNCTION TO CREATE A SPATIAL GRID
# 'xx' have a 'lat' and a 'long' columns
assign.points.to.a.spatial.grid <- function(xx, general){

  xx <- xx[,!colnames(xx) %in% c("icessquare","icessquare.vms") ]  # remove
  xx <- cbind.data.frame(xx, icessquare= rep(0,nrow(xx)))


  an <- function(x) as.numeric(as.character(x))
  rlong      <- range(an(xx$long),na.rm=T)
  vect.long  <- signif(seq(floor(rlong[1]), ceiling(rlong[2]), by=1),4)   # long (x)
  label.long <- rep(paste(rep(LETTERS,each=10),0:9,sep=""),each=1)
  names(label.long) <- signif(seq(-50, 209, by=1),4)   # long (x)
  label.long <- label.long[!is.na(names(label.long))]  # => correspondance long (-50 to 209) / sq letter (A0 to Z9)
  label.long <- label.long[as.character(vect.long)]
  rlat      <- range(an(xx$lat), na.rm=T)
  vect.lat   <- signif(seq(floor(rlat[1]), ceiling(rlat[2]),by=0.5),4) # lat  (y)
  label.lat  <- rep(paste(seq(1,75,1)),each=1)
  names(label.lat) <-   paste(signif(seq(36,73, by=0.5),4))
  label.lat <- label.lat[!is.na(names(label.lat))] # => correspondance lat (36 to 73) / sq number (1 to 75)
  label.lat <- label.lat[as.character(vect.lat)]
  vect.label <- paste(rep(label.lat,each=length(label.long)),"",label.long,sep="")
  xx[,"icessquare"] <- paste(label.lat [findInterval(an(xx[,"lat"]), vect.lat)] , label.long [findInterval(an(xx[,"long"]), vect.long)], sep="")

  return(xx)
  }


#!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!#
#utils--
# for managing NA on logbook side
# (from vms trip.sq without corresponding logbook trip.sq e.g. because no declaration in sq because only steaming time inside)
# we need to inform back the specificity of the vessel from logbook using info from the same trip i.e. vesselid+bk.tripnum
 retrieve.on.bk.side <- function(merged, type.data){
      idx <- which(merged$metier=="NA")
      merged.NA <- merged[idx,] # input (only the trip.sq with NA for the logbook part)

      for (td in type.data){
         map <- tapply(merged[, td ], paste(merged$vesselid, merged$bk.tripnum),
                             function(i) {ss<- unique(as.character(i)) ; ss[ss!="NA"][1]})
         merged.NA[, td ] <- factor(paste(merged.NA$vesselid,merged.NA$bk.tripnum))
         levels(merged.NA[, td ]) <- map[levels(merged.NA[, td ])]
         }
      if(nrow(merged.NA)>0) merged.NA$flag <- 4 # flag on meth
      merged[idx,] <- merged.NA # output
      return(merged)
      }



#!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!#
#utils--
  export.check.merging.quality <- function(a.case=1,  general=general, 
                      print.console=TRUE, is.check=TRUE,  fuelcons=FALSE, ...){
                arg <- list(...)
         if(a.case==1){
                if(is.check){
                  # vesselid
                  a.vesselid <- as.character(arg$agg.logbk.this.vessel.method.3[1,"vesselid"])
                  # year
                  a.year <- general$a.year
                  # method
                  a.method <- general$landings.redistribution
                  # weight sp.to.keep
                  nm      <- names(arg$agg.logbk.this.vessel.method.3)
                  idx.col.w <- grep('KG', nm) # index columns with species weight
                  idx.col.v <- grep('EURO', nm) # index columns with species value
                  idx.col <- c(idx.col.w, idx.col.v)
                  weight.per.sp <- apply(arg$agg.logbk.this.vessel.method.3[,idx.col], 2, sum, na.rm=TRUE)
                  weight.per.sp.to.keep <- weight.per.sp[general$sp.to.keep]
                }
                inp.w <- sum(arg$agg.logbk.this.vessel.method.3[, nm %in% general$sp.to.keep[1] ])
                if(print.console){
                   cat(paste("input weight in logbk for", general$sp.to.keep[1],":", signif(inp.w,5), "\n",sep=' '))
                   assign("inp.w", inp.w, envir=sys.frame(-1))
                   }
                # effort state1, state2
                if(general$landings.redistribution=="mixture123"){
                        inp.v <- sum(arg$.vms$effort.mins)
                        inp.f <- sum(arg$.vms$fuelcons)
                        if(print.console) cat(paste("input effort in vms:", signif(inp.v,5), "\n",sep=' '))
                        if(fuelcons) if(print.console) cat(paste("input fuelcons in vms:", signif(inp.f,5), "\n",sep=' '))
                if(is.check){
                  input.effort.state1 <-  sum(arg$.vms[arg$.vms$state==1,]$effort.mins)
                  input.effort.state2 <-  sum(arg$.vms[arg$.vms$state==2,]$effort.mins)
                  }
                }
                if(print.console) assign("inp.v", inp.v, envir=sys.frame(-1))
                if(fuelcons) if(print.console) assign("inp.f", inp.f, envir=sys.frame(-1))
                # fill in the 'exported.vector'
                if(is.check) exported.vector [c("vesselid","year", "method",
                      paste("input.weight.",general$sp.to.keep,sep=''), "input.effort.state1", "input.effort.state2")]  <<- c(
                        a.vesselid, a.year, a.method, weight.per.sp.to.keep, input.effort.state1, input.effort.state2)
           }
           if(a.case==2){
                      if(is.check) {
                         nm      <- names(arg$.logbk.for.meth1)
                         idx.col.w <- grep('KG', nm) # index columns with species weight
                         idx.col.v <- grep('EURO', nm) # index columns with species value
                         idx.col <- c(idx.col.w, idx.col.v)
                         weight.per.sp <- apply(arg$.logbk.for.meth1[,idx.col], 2, sum, na.rm=TRUE)
                         weight.per.sp.to.keep <- weight.per.sp[general$sp.to.keep]
                         }
                     if(print.console) {
                         inp <- sum(arg$.logbk.for.meth1[, nm %in% general$sp.to.keep[1] ])
                         cat(paste("input weight in logbk for ",general$sp.to.keep[1]," for meth1:", signif(inp,5), "\n",sep=' '))
                         }
                      if(general$landings.redistribution=="mixture123"){
                          if(is.check) {
                           input.effort.state1 <-  sum(arg$.vms.for.meth1[arg$.vms.for.meth1$state==1,]$effort.mins)
                           input.effort.state2 <-  sum(arg$.vms.for.meth1[arg$.vms.for.meth1$state==2,]$effort.mins)
                            }
                          if(print.console) cat(paste("input effort in vms for meth1:", signif(sum(arg$.vms.for.meth1$effort.mins),5), "\n",sep=' '))
                       }
                      if(general$landings.redistribution=="mixture456"){
                         dd <- arg$.vms.for.meth1[!duplicated(data.frame(arg$.vms.for.meth1$icessquare.vms,arg$.vms.for.meth1$bk.tripnum)),]
                         if(is.check){
                            input.effort.state1 <- sum(an(dd$effort.state1), na.rm=T)
                            input.effort.state2 <-  sum(an(dd$effort.state2), na.rm=T)
                            }
                         if(print.console) {
                            inp <- sum(an(dd$effort.state2), na.rm=T)+ sum(an(dd$effort.state1), na.rm=T) # => input
                            cat(paste("input effort in vms:", signif(inp,5), "\n",sep=' '))
                            }
                         }
                # fill in the 'exported.vector'
                if(is.check) exported.vector [c(paste("input.weight.meth1.",general$sp.to.keep,sep=''),
                                     "input.effort.state1.meth1", "input.effort.state2.meth1")]  <<- c(
                                weight.per.sp.to.keep, input.effort.state1, input.effort.state2)
           }
           if(a.case==3){
                          if(is.check){
                          nm      <- names(arg$.logbk.for.meth2)
                         idx.col.w <- grep('KG', nm) # index columns with species weight
                         idx.col.v <- grep('EURO', nm) # index columns with species value
                         idx.col <- c(idx.col.w, idx.col.v)
                         weight.per.sp <- apply(arg$.logbk.for.meth2[,idx.col], 2, sum, na.rm=TRUE)
                         weight.per.sp.to.keep <- weight.per.sp[general$sp.to.keep]
                         }
                         if(print.console){
                         inp <- sum(arg$.logbk.for.meth2[, nm %in% general$sp.to.keep[1] ])
                          cat(paste("input weight in logbk for", general$sp.to.keep[1], "for meth2:", signif(inp,5),"\n", sep=' '))
                          }
                         if(general$landings.redistribution=="mixture123"){
                            if(print.console) cat(paste("input effort in vms for meth2:", signif(sum(arg$.vms.for.meth2$effort.mins),5), "\n",sep=' '))
                            if(is.check){
                              input.effort.state1 <-  sum(arg$.vms.for.meth2[arg$.vms.for.meth2$state==1,]$effort.mins)
                              input.effort.state2 <-  sum(arg$.vms.for.meth2[arg$.vms.for.meth2$state==2,]$effort.mins)
                            }

                         }
                         if(general$landings.redistribution=="mixture456"){
                           dd <- arg$.vms.for.meth2[!duplicated(data.frame(arg$.vms.for.meth2$icessquare.vms,arg$.vms.for.meth2$bk.tripnum)),]
                           if(print.console) {
                                inp <- sum(an(dd$effort.state2), na.rm=T)+ sum(an(dd$effort.state1), na.rm=T) # => input
                                cat(paste("input effort in vms:", signif(inp,5), "\n",sep=' '))
                           }
                           if(is.check){
                                input.effort.state1 <- sum(an(dd$effort.state1), na.rm=T)
                                input.effort.state2 <-  sum(an(dd$effort.state2), na.rm=T)
                           }
                         }
                # fill in the 'exported.vector'
                if(is.check) exported.vector [c(paste("input.weight.meth2.",general$sp.to.keep,sep=''),
                                     "input.effort.state1.meth2", "input.effort.state2.meth2")]  <<- c(
                                       weight.per.sp.to.keep, input.effort.state1, input.effort.state2)
               }
           if(a.case==4){
                          if(is.check){
                          nm      <- names(arg$.logbk.for.meth3)
                         idx.col.w <- grep('KG', nm) # index columns with species weight
                         idx.col.v <- grep('EURO', nm) # index columns with species value
                         idx.col <- c(idx.col.w, idx.col.v)
                         weight.per.sp <- apply(arg$.logbk.for.meth3[,idx.col], 2, sum, na.rm=TRUE)
                         weight.per.sp.to.keep <- weight.per.sp[general$sp.to.keep]
                         }
                         if(print.console){
                         inp <- sum(arg$.logbk.for.meth3[, nm %in% general$sp.to.keep[1] ])
                          cat(paste("input weight in logbk for", general$sp.to.keep[1], "for meth3:", signif(inp,5),"\n", sep=' '))
                          }
                         if(general$landings.redistribution=="mixture123"){
                            if(print.console) cat(paste("input effort in vms for meth3:", signif(sum(arg$.vms.for.meth3$effort.mins),5), "\n",sep=' '))
                            if(is.check){
                              input.effort.state1 <-  sum(arg$.vms.for.meth3[arg$.vms.for.meth3$state==1,]$effort.mins)
                              input.effort.state2 <-  sum(arg$.vms.for.meth3[arg$.vms.for.meth3$state==2,]$effort.mins)
                            }

                         }
                 # fill in the 'exported.vector'
                if(is.check) exported.vector [c(paste("input.weight.meth3.",general$sp.to.keep,sep=''),
                                     "input.effort.state1.meth3", "input.effort.state2.meth3")]  <<- c(
                                       weight.per.sp.to.keep, input.effort.state1, input.effort.state2)
               }
           if(a.case==5){
                   idx <- which(arg$merged$metier=="NA")
                   merged.NA <- arg$merged[idx,]
                   dd <- merged.NA[!duplicated(data.frame(merged.NA$idx)),]
                   if(general$landings.redistribution=="mixture123"){
                       if(print.console) {
                                  out.v <- sum(an(dd$effort.mins), na.rm=TRUE)
                                  cat(paste("vms effort not in bk:", signif(out.v,5), "\n", sep=' '))
                        }
                       if(is.check){
                        effort.not.bk.state1 <- sum(an(dd[dd$state==1,]$effort.mins), na.rm=TRUE)
                        effort.not.bk.state2 <- sum(an(dd[dd$state==2,]$effort.mins), na.rm=TRUE)
                      }
                   }
               # fill in the 'exported.vector'
                if(is.check) exported.vector [c("effort.not.bk.state1", "effort.not.bk.state2")]  <<- c(
                        effort.not.bk.state1, effort.not.bk.state2)

           }
           if(a.case==6){
                      if(is.check){
                         nm      <- names(arg$merged)
                         idx.col.w <- grep('KG', nm) # index columns with species weight
                         idx.col.v <- grep('EURO', nm) # index columns with species value
                         idx.col <- c(idx.col.w, idx.col.v)
                         weight.per.sp <- apply(arg$merged[,idx.col], 2, sum, na.rm=TRUE)
                         weight.per.sp.to.keep <- weight.per.sp[general$sp.to.keep]
                         }
                      if(print.console){
                        out.w <- sum(arg$merged[, nm %in% general$sp.to.keep[1] ], na.rm=TRUE)
                         cat(paste("output weight in 'merged' for", general$sp.to.keep[1],":", signif(out.w,5), sep=' '))
                         if(out.w > arg$inp.w*1.1 || out.w < arg$inp.w*0.9)  cat(" <<\n") else cat("\n")
                         }
                      if(general$landings.redistribution=="mixture123"){
                              dd <- arg$merged[!duplicated(data.frame(arg$merged$idx)),]
                              if(print.console) {
                                  out.v <- sum(an(dd$effort.mins), na.rm=TRUE)
                                  cat(paste("output effort in 'merged':", signif(out.v,5), sep=' '))
                                  if(fuelcons) out.f <- sum(an(dd$fuelcons), na.rm=TRUE)
                                  if(out.v > arg$inp.v*1.1 || out.v < arg$inp.v*0.9)  cat(" <<\n") else cat("\n")
                                  if(fuelcons) cat(paste("output fuelcons in 'merged':", signif(out.f,5), sep=' '))
                                  if(fuelcons) if(out.f > arg$inp.f*1.1 || out.f < arg$inp.f*0.9) cat(" <<\n") else  cat("\n")
                               }
                             if(is.check){
                              output.effort.state1 <-  sum(an(dd[dd$state==1,]$effort.mins), na.rm=TRUE)
                              output.effort.state2 <-  sum(an(dd[dd$state==2,]$effort.mins), na.rm=TRUE)
                              }
                         }
                   # fill in the 'exported.vector'
                if(is.check) exported.vector [c(paste("output.weight.merged.",general$sp.to.keep,sep=''),
                                     "output.effort.state1", "output.effort.state2")]  <<- c(
                        weight.per.sp.to.keep, output.effort.state1, output.effort.state2)
           }
           if(a.case=="export.this.vessel"){
                   # export
                  if(is.check) write.table(t(exported.vector), file= file.path(general$main.path,
                     paste("export-check-merging-quality-eflalo-",general$a.year,".txt",sep="")),
                       append=TRUE, row.names=FALSE, col.names=FALSE)
                   }


 return()
 }
 
 
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
merge.vms.to.logbook.at.the.ping.scale <-
           function(logbooks, vms, general=general, ...){

an <<- function(x) as.numeric(as.character(x)) # alias

  lstargs <- list(...)
  an <- function(x) as.numeric(as.character(x))

  # (optional)
  if(!is.null(general$export.check.merging.quality) && general$export.check.merging.quality){
     general$nm.exported.vector <-
             c("vesselid","year", "method",
               paste("input.weight.",general$sp.to.keep,sep=''), "input.effort.state1", "input.effort.state2",
                paste("input.weight.meth1.",general$sp.to.keep,sep=''), "input.effort.state1.meth1", "input.effort.state2.meth1",
                  paste("input.weight.meth2.",general$sp.to.keep,sep=''), "input.effort.state1.meth2", "input.effort.state2.meth2",
                    paste("input.weight.meth3.",general$sp.to.keep,sep=''), "input.effort.state1.meth3", "input.effort.state2.meth3",
                   "effort.not.bk.state1", "effort.not.bk.state2",
                    paste("output.weight.merged.",general$sp.to.keep,sep=''), "output.effort.state1", "output.effort.state2")
     write.table(t(general$nm.exported.vector), file=file.path(general$main.path,  # init
                   paste("export-check-merging-quality-eflalo-",general$a.year,".txt",sep="")), append=FALSE, row.names=FALSE, col.names=FALSE)
     }
     # for exporting the merging quality check for this vessel
     exported.vector <<- NULL


      #!#!##!#!##!#!##!#!##!#!##!#!#
      #!#!##!#!##!#!##!#!##!#!##!#!#
      #!#!##!#!##!#!##!#!##!#!##!#!#
      #!#!##!#!##!#!##!#!##!#!##!#!#
      #!#!##!#!##!#!##!#!##!#!##!#!#
      all.vesselid     <- as.character(unique(logbooks[an(logbooks$VE_LEN)>1,]$VE_REF)) 
      # => keep vessels > 15 m (i.e. vms-equipped)
      if(length(lstargs$a.vesselid)!=0) 
          all.vesselid <- lstargs$a.vesselid # KEEP ONLY ONE VESSEL IF NEEDED....

      for(a.vesselid in all.vesselid){  # PER VESSEL

         #----------
         #----------
         #----------
         #----------
         #----------
         #----------
         # LOGBOOK INPUT
         logbk.this.vessel            <- logbooks[logbooks$VE_REF %in% a.vesselid,]
         logbk.this.vessel$LE_RECT    <- factor(logbk.this.vessel$LE_RECT)
         logbk.this.vessel$VE_REF     <- factor( logbk.this.vessel$VE_REF)

         #  if it still does not exist, add mid-time and bk.tripnum in logbook
        if(!any(colnames(logbk.this.vessel)%in%c('mid.time','bk.tripnum'))){
           # departure time
           Sys.setlocale("LC_TIME", "english")
           ctime <- strptime(  paste(logbk.this.vessel$FT_DDAT, logbk.this.vessel$FT_DTIME) ,  "%e/%m/%Y %H:%M" )
           logbk.this.vessel <- cbind.data.frame(logbk.this.vessel, date.in.R.dep=ctime)
           # arrival time
           Sys.setlocale("LC_TIME", "english")
           ctime <- strptime(  paste(logbk.this.vessel$FT_LDAT, logbk.this.vessel$FT_LTIME) ,  "%e/%m/%Y %H:%M" )
           logbk.this.vessel <- cbind.data.frame(logbk.this.vessel, date.in.R.arr=ctime)
           # catch.date
           Sys.setlocale("LC_TIME", "english")
           ctime <- strptime(  paste(logbk.this.vessel$LE_CDAT) ,  "%e/%m/%Y" )
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
           }

         if(!is.null(general$metier.def)){
           # metier definition 
           logbk.this.vessel$LE_MET   <- factor( logbk.this.vessel$VE_REF)
           if(general$metier.def=="gear_meshsize") {
             logbk.this.vessel$LE_MET <-
               paste(logbk.this.vessel$LE_GEAR, logbk.this.vessel$LE_MSZ,sep="_")
           }
           if(general$metier.def=="gear_meshsize_targetpca") {
             logbk.this.vessel$LE_MET <-
                paste(logbk.this.vessel$LE_GEAR, logbk.this.vessel$LE_MSZ, logbk.this.vessel$target.trip.pca, sep="_")
           }
           if(general$metier.def=="gear_targetpca") {
             logbk.this.vessel$LE_MET <-
                paste(logbk.this.vessel$LE_GEAR, logbk.this.vessel$target.trip.pca, sep="_")
           }
         }
         #=> LOGBOOK INPUT REQUIRES AT LEAST,
         #     'VE_REF',  FT_DDAT, FT_DTIME, FT_LDAT, FT_LTIME, FT_CDAT,
         #  'LE_SP_KG' (etc.), 'LE_RECT', 'VE_FLT' AND 'LE_MET' (informed or not) AND 'LE_GEAR' COLUMNS
         #

         #----------
         #----------
         #----------
         #----------
         #----------
         #----------
         # VMS INPUT: load traj with 'at sea' pings
         # ABSOLUTELY REQUIRED: c("VE_REF","SI_LATI","SI_LONG", "SI_DATE", "SI_TIME", "SI_FT", "SI_HARB", "SI_STATE")
         # optional: "fuelcons"


         if(a.vesselid %in% unique(tacsat$VE_REF)){
         
         tacsat <- tacsat[tacsat$VE_REF == a.vesselid,] # subset for this vessel
         tacsat$VE_REF <- factor(tacsat$VE_REF)
         
         # convert TACSAT names in local names (if required)
         colnames(tacsat)  [colnames(tacsat) %in% "VE_REF"]   <- "vesselid"
         colnames(tacsat)  [colnames(tacsat) %in% "SI_FT"]    <- "tripnum"
         colnames(tacsat)  [colnames(tacsat) %in% "SI_LATI"]  <- "lat"
         colnames(tacsat)  [colnames(tacsat) %in% "SI_LONG"]  <- "long"
         colnames(tacsat)  [colnames(tacsat) %in% "SI_STATE"] <- "state"
         colnames(tacsat)  [colnames(tacsat) %in% "SI_HARB"]  <- "which"

         # if does not exist, add date.in.R for handling time
         if(!("date.in.R" %in% colnames(tacsat))){
           Sys.setlocale("LC_TIME", "english")
           ctime <- strptime(  paste(tacsat$SI_DATE, tacsat$SI_TIME) ,  "%e/%m/%Y %H:%M" )
           tacsat <- cbind.data.frame(tacsat, date.in.R=ctime)
         }

         # keep only the essential
         if("fuelcons" %in% colnames(tacsat)){
            vms.this.vessel  <- tacsat [, c("vesselid","lat","long", "date.in.R","tripnum", "which", "state", "fuelcons")]
         }else{
            vms.this.vessel  <- tacsat [, c("vesselid","lat","long", "date.in.R","tripnum", "which", "state")]
         }


         vms.this.vessel$idx  <- 1:nrow(vms.this.vessel) # label for each ping

         vms.this.vessel$vesselid   <- factor(vms.this.vessel$vesselid)

         # filter if vessel with a bad vms
         to.remove.because.deficient.vms <- any(is.na(vms.this.vessel$tripnum))
         to.remove.because.not.enough.vms.trips <- length(table(vms.this.vessel$tripnum))< 5  # nb vms trips < 5
         flag <- to.remove.because.deficient.vms ||  to.remove.because.not.enough.vms.trips

         ## remove bk.tripnum and mid.time if it exists
         vms.this.vessel <- vms.this.vessel[, !colnames(vms.this.vessel) %in% c("bk.tripnum", "mid.time")]



        if(flag==FALSE) {  # i.e. vms-equipped


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
           .vms$start.trip <- c(1,diff(.vms[,"tripnum"]))
           .vms$end.trip <- c(diff(.vms[,"tripnum"]),0)
           .vms[.vms$start.trip>0, "start.trip"] <- .vms[.vms$start.trip>0, "tripnum"]
           .vms[.vms$end.trip>0, "end.trip"] <- .vms[.vms$end.trip>0, "tripnum"]

           tmp <- .vms[.vms$start.trip>0,]
           tmp <- tmp[,c("vesselid","date.in.R","tripnum")]
           tmp2 <- .vms[.vms$end.trip>0,]
           tmp2 <- tmp2[,c("vesselid","date.in.R","tripnum")]
           .vms <- .vms[,!colnames(.vms) %in% c("start.trip", "end.trip")] # remove tool columns
           table.midtime <- merge(tmp, tmp2, by.x="tripnum", by.y="tripnum")
           table.midtime <- table.midtime[, c("tripnum","vesselid.x","date.in.R.x","date.in.R.y") ]
           colnames(table.midtime) <- c("tripnum","vesselid","date.in.R.dep","date.in.R.arr")
           } else{
           table.midtime <- .vms[, c("tripnum","vesselid","date.in.R.dep","date.in.R.arr") ]
           table.midtime <- table.midtime[!duplicated(data.frame(table.midtime$tripnum, table.midtime$vesselid)),]
           }
         } else{stop("no 'date.in.R' found in vms")}
         mid.time <- rep(0, nrow(table.midtime))
         for(r in 1: nrow(table.midtime)){
           mid.time[r] <-  as.character(seq(from=table.midtime$date.in.R.dep[r], to=table.midtime$date.in.R.arr[r], length.out = 3)[2])
         }
         table.midtime$mid.time <-  mid.time
         if(!any(colnames(.vms)%in%"mid.time")){ # here we are...
              .vms <- merge(.vms, table.midtime[,c("tripnum","mid.time")], by.x="tripnum", by.y="tripnum")
         }

        #!!!!!!!!!!!!!!!!!!#
        #!!!!!!!!!!!!!!!!!!#
        # -IF IT DOES NOT EXIST YET-,
        # ASSIGN A 'BK.TRIPNUM' FROM LOGBOOK TO EACH VMS TRIP
       if(!any(colnames(.vms)%in%"bk.tripnum")){
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
              text(as.POSIXct(table.midtime$mid.time[i]), 0.52, table.midtime$tripnum[i], cex=0.5, col=1)
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
          # (so, for each mid.time in vms, a tripnum will be find)
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
            sov <- as.character(.logbk$VE_FLT[1])
            savePlot(filename = file.path(general$main.path,
                            paste("assign_bk.tripnum_to_vms_",ve,"-",sov,"-",general$a.year,".jpeg",sep="")),type ="jpeg")
           dev.off()
          }

       } else{
         if(TRUE) { # if already done, just remenber the first merging
            the.vms.ve <- .vms[!duplicated(data.frame(.vms$bk.tripnum, .vms$vesselid)),]
            the.logbk.ve <- .logbk[!duplicated(data.frame(.logbk$bk.tripnum, .logbk$VE_REF)),]

            trunk <- 1 # trunk give the part of the year to be plotted (1 to 5)
            windows(width=8, height=4)
            ltrunk <- (nrow(the.vms.ve)/3)
            idxtrunk <-  (trunk+(trunk-1)*ltrunk):(trunk*ltrunk)

            plot(the.vms.ve$date.in.R.dep[ idxtrunk],rep(1,length(the.vms.ve$date.in.R.dep[ idxtrunk])),
                 ylim=c(0,0.52), type="n", ylab="", axes=FALSE)
            r <- as.POSIXct(round(range(the.vms.ve$date.in.R.dep), "days"))
            axis.POSIXct(1, at=seq(r[1], r[2], by="month"), format="%e%b%y:%H:%M")
            axis(2, at=c(0.5,0.1),labels=c("VMS","LOGBOOK"))
            for(i in 1:nrow(the.vms.ve))  {
               segments(as.POSIXct(the.vms.ve$date.in.R.dep[i]), 0.5, as.POSIXct(the.vms.ve$date.in.R.arr[i]), 0.5, col=1)
               tmp <- the.logbk.ve[the.logbk.ve$bk.tripnum  %in% the.vms.ve$bk.tripnum[i],]$mid.time
               #if(length(tmp!=0))
               arrows(as.POSIXct(the.vms.ve$date.in.R.dep[i] ), 0.5 ,as.POSIXct(tmp ),0.1, length=0.1)
            }
            for(i in 1:nrow(the.logbk.ve))  {
               segments(as.POSIXct(the.logbk.ve$date.in.R.dep[i]), 0.1, as.POSIXct(the.logbk.ve$date.in.R.arr[i]), 0.1, col=1)
               text(as.POSIXct(the.logbk.ve$date.in.R.dep[i]), jitter(0.05,amount=0.01), as.character(the.logbk.ve$bk.tripnum[i]), cex=0.5)
            }
     } # end FALSE
    } # end else
         .logbk$mid.time    <- factor(.logbk$mid.time)
         .logbk$bk.tripnum  <- factor(.logbk$bk.tripnum)
         .vms$mid.time      <- factor(.vms$mid.time)
         .vms$bk.tripnum    <- factor(.vms$bk.tripnum)

         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         #! ASSIGN A 'TRIPNUM' FROM VMS TRIP NUM TO  #!#!#
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
         .vms   <- assign.points.to.a.spatial.grid(xx=.vms, general)
        
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         # COMPUTE EFFORT.MINS      !#!#!#!#!#!#!#!#!#!#!#
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
          .vms <- .vms[order(.vms$date.in.R),]
          .vms$effort.mins <- abs(c(0, as.numeric(.vms[-nrow(.vms),"date.in.R"] - 
                                        .vms[-1,"date.in.R"], units="mins")))
           start.trip <- c(1,diff(.vms[,"tripnum"]))
          .vms[start.trip!=0, "effort.mins"] <- 0  # just correct for the trip change points



         rm(er); rm(xx) ; gc(reset=TRUE)
                                           
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         # SET UP PRIMARY KEYS FOR MERGING!#!#!#!#!#!#!#!#
         #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
         .logbk$bk.tripnum <- factor(.logbk$bk.tripnum )
         .logbk$bk.tripnum.sq <- paste(.logbk$bk.tripnum, ".", .logbk$LE_RECT, sep='') # caution:redefine
         .logbk$bk.tripnum.sq.day <- paste(.logbk$bk.tripnum, ".", .logbk$LE_RECT,".",.logbk$date.in.R.cat, sep='') # caution:redefine
         .vms$bk.tripnum <- factor(.vms$bk.tripnum)
         .vms$bk.tripnum.sq <- paste(.vms$bk.tripnum, ".", .vms$icessquare, sep='') # caution:redefine
         .vms$bk.tripnum.sq.day <- paste(.vms$bk.tripnum, ".", .vms$icessquare,".", format(.vms$date.in.R,  '%Y-%m-%d'), sep='') # caution:redefine

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
                              .logbk$VE_REF, .logbk$VE_FLT,  .logbk$LE_MET, .logbk$LE_GEAR), sum, na.rm=TRUE )
              colnames(agg.logbk.this.vessel.method.1) <- 
                           c("bk.tripnum", "vesselid", "fleet", "metier","gear", nm[idx.col] )
             # AGGREGATE WEIGHT (OR VALUE) PER SPECIES PER BK.TRIPNUM.SQ
              agg.logbk.this.vessel.method.2  <- aggregate(.logbk[,idx.col],
                      list(.logbk$bk.tripnum.sq, 
                              .logbk$VE_REF, .logbk$VE_FLT,  .logbk$LE_MET, .logbk$LE_GEAR), sum, na.rm=TRUE )
              colnames(agg.logbk.this.vessel.method.2) <- 
                           c("bk.tripnum.sq", "vesselid", "fleet","metier" ,"gear", nm[idx.col])
             # AGGREGATE WEIGHT (OR VALUE) PER SPECIES PER BK.TRIPNUM.SQ.DAY (NOTE: SO, 'LE_SEQNUM' IS AGGREGATED HERE)
              agg.logbk.this.vessel.method.3  <- aggregate(.logbk[,idx.col],
                      list(.logbk$bk.tripnum.sq.day, 
                             .logbk$VE_REF, .logbk$VE_FLT,  .logbk$LE_MET, .logbk$LE_GEAR), sum, na.rm=TRUE )
              colnames(agg.logbk.this.vessel.method.3) <- 
                          c("bk.tripnum.sq.day", "vesselid", "fleet","metier","gear",  nm[idx.col])


             #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
             # MERGING WITH VMS PER TRIP !!!!!!!!!!#!#!#!#!#!#
             #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
             do.merging <- function(method="bk.tripnum", .logbk, .vms, general){


              # IF BY PING-------------
              # find total nb of FISHING ping per tripnum from vms  # used for method 1  'bk.tripnum'
              if(method=="bk.tripnum" &&
                   (general$landings.redistribution=="method1" || general$landings.redistribution=="method2" ||
                      general$landings.redistribution=="method3" ||
                       general$landings.redistribution=="mixture123")){
               .vms$count.fping.trip  <- factor(.vms$bk.tripnum)  # init
              count.fping.trip <- table(.vms[.vms$state==1,]$bk.tripnum)
              # => COUNT nb of FISHING pings per bk.tripnum because each weight will be repeated by ping after merging
              levels(.vms$count.fping.trip) <- count.fping.trip[levels(.vms$count.fping.trip)]  # mapping
              .vms[.vms$state==2,]$count.fping.trip <- NA
              # => COUNT nb of gears per bk.tripnum because each ping will be repeated by gear after merging
              count.gr.trip <- tapply(.logbk$gear, .logbk$bk.tripnum, function(x) length(unique(x)))
              .logbk$count.gr.trip <- count.gr.trip[.logbk$bk.tripnum]  # mapping
              }


              # find total nb of FISHING ping per trip-icessquare from vms  # used for method 2   'bk.tripnum.sq'
              if(method=="bk.tripnum.sq" &&
                   (general$landings.redistribution=="method1" || general$landings.redistribution=="method2" ||
                      general$landings.redistribution=="method3" ||
                       general$landings.redistribution=="mixture123")){
              .vms$count.fping.trip.sq  <- factor(.vms$bk.tripnum.sq)  # init
              count.fping.trip.sq <- table(.vms[.vms$state==1,]$bk.tripnum.sq) # COUNT nb of FISHING pings per bk.tripnum.sq
              levels(.vms$count.fping.trip.sq) <- count.fping.trip.sq[levels(.vms$count.fping.trip.sq)]  # mapping
              if(any('2' %in% unique(.vms$state))) .vms[.vms$state==2,]$count.fping.trip.sq <- NA
              # => COUNT nb of gears per bk.tripnum.sq because each ping will be repeated by gear after merging
              count.gr.trip.sq <- tapply(.logbk$gear, .logbk$bk.tripnum.sq, function(x) length(unique(x)))
              .logbk$count.gr.trip.sq <- count.gr.trip.sq[.logbk$bk.tripnum.sq]  # mapping
              }

              # find total nb of FISHING ping per trip-icessquare-day from vms  # used for method 3   'bk.tripnum.sq.day'
              if(method=="bk.tripnum.sq.day" &&
                   (general$landings.redistribution=="method1" || general$landings.redistribution=="method2" ||
                      general$landings.redistribution=="method3" ||
                       general$landings.redistribution=="mixture123")){
              .vms$count.fping.trip.sq.day  <- factor(.vms$bk.tripnum.sq.day)  # init
              count.fping.trip.sq.day <- table(.vms[.vms$state==1,]$bk.tripnum.sq.day) # COUNT nb of FISHING pings per bk.tripnum.sq.day
              levels(.vms$count.fping.trip.sq.day) <- count.fping.trip.sq.day[levels(.vms$count.fping.trip.sq.day)]  # mapping
              if(any('2' %in% unique(.vms$state))) .vms[.vms$state==2,]$count.fping.trip.sq.day <- NA
              # => COUNT nb of gears per bk.tripnum.sq.day because each ping will be repeated by gear after merging
              count.gr.trip.sq.day <- tapply(.logbk$gear, .logbk$bk.tripnum.sq.day, function(x) length(unique(x)))
              .logbk$count.gr.trip.sq.day <- count.gr.trip.sq.day[.logbk$bk.tripnum.sq.day]  # mapping}
              }



              # do the merging between .logbk and .vms according to
              #  meth1: 'bk.tripnum' OR meth2: 'bk.tripnum.sq' OR meth3: 'bk.tripnum.sq.day'
              # need to use a trick to avoid "out of memory" doing the merge()
              names(.logbk)  [names(.logbk) %in% "VE_REF"] <- "vesselid" # change needed to match with .vms
              coln.idx1 <- which(!colnames(.logbk)%in%c("vesselid", method))
              coln1 <- colnames(.logbk)[coln.idx1]
              tmp1 <- data.frame(coll= collapse.all.columns  (.logbk, columns= coln.idx1  ),
                         vesselid=.logbk$vesselid, a.method= .logbk[,method] ) #.logbk
              coln.idx2 <- which(!colnames(.vms)%in%c("vesselid", method))
              coln2 <- colnames(.vms)[coln.idx2]
              tmp2 <- data.frame(coll2= collapse.all.columns  (.vms, columns=  coln.idx2 ),
                         vesselid=.vms$vesselid, a.method= .vms[,method] )  #.vms
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
              tmp5 <- cbind.data.frame(merged.this.vessel[,c("vesselid", method)], tmp3, tmp4)
              colnames(tmp5) <- c("vesselid", method, coln1, coln2)
              merged.this.vessel <- tmp5

              # note: at this stage some few lgbk records could have got 
              # no correpondance in vms (misreporting of area)=> NA on vms side
              # and then the landing weight of these records are possibly lost because count.ping at NA...
              # so we can choose to correct (see ** below) to keep these land. weight
              # the remaining loss in weight will come from the matching records having catches but
              # without fishing pings (i.e. only steaming pings)!

          

              if(TRUE){
              # conservation of catches?
              # detect possible weight landed while no feffort detected from vms
              if(general$landings.redistribution=="mixture123"){
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
              }
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
              if(method=="bk.tripnum.sq.day" && (general$landings.redistribution=="method3" ||
                                   general$landings.redistribution=="mixture123") ){
                  merged.this.vessel[merged.this.vessel$lat=='NA', "count.fping.trip.sq.day"] <- 1 # **correct for loss if in lgk but not in vms
                             merged.this.vessel[,idx.col] <- (apply(merged.this.vessel[,idx.col],2,an) /
                                                        an(merged.this.vessel$count.fping.trip.sq.day)) /
                                                                        an(merged.this.vessel$count.gr.trip.sq.day)
              }
              if(method=="bk.tripnum.sq" && (general$landings.redistribution=="method2" ||
                                   general$landings.redistribution=="mixture123") ){
                  merged.this.vessel[merged.this.vessel$lat=='NA', "count.fping.trip.sq"] <- 1 # **correct for loss if in lgk but not in vms
                             merged.this.vessel[,idx.col] <- (apply(merged.this.vessel[,idx.col],2,an) /
                                                        an(merged.this.vessel$count.fping.trip.sq)) /
                                                                        an(merged.this.vessel$count.gr.trip.sq)
              }
              if(method=="bk.tripnum" && (general$landings.redistribution=="method1" ||
                                   general$landings.redistribution=="mixture123") ){
                 merged.this.vessel[merged.this.vessel$lat=='NA', "count.fping.trip"] <- 1 # **correct for loss if in lgk but not in vms
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
             if(general$export.check.merging.quality) is.check <- TRUE else is.check <- FALSE
              if(general$landings.redistribution=="method1"){
              .logbk   <- agg.logbk.this.vessel.method.1 # METHOD 1: per tripnum
                merged <- do.merging(method="bk.tripnum", .logbk, .vms, general)
              }
              if(general$landings.redistribution=="method2" ){
              .logbk   <- agg.logbk.this.vessel.method.2 # METHOD 2: per tripnum.sq
                merged <- do.merging(method="bk.tripnum.sq", .logbk, .vms, general)
              }
              if(general$landings.redistribution=="method3"){
              .logbk   <- agg.logbk.this.vessel.method.3 # METHOD 3: per tripnum.sq.day
                merged <- do.merging(method="bk.tripnum.sq.day", .logbk, .vms, general)
              }
              if(general$landings.redistribution=="mixture123"){
                 .logbk   <- agg.logbk.this.vessel.method.3
                 my.split <- function(obj,a.sep="\\.",idx=1) unlist(lapply(strsplit(obj, a.sep),function(x)x[idx]))
                 # reduce the level
                 .logbk$bk.tripnum.sq  <-  paste(my.split(as.character(.logbk$bk.tripnum.sq.day),a.sep="\\.",idx=1),
                                                 my.split(as.character(.logbk$bk.tripnum.sq.day),a.sep="\\.",idx=2),sep='.')
                 # reduce the level
                 .logbk$bk.tripnum     <-        my.split(as.character(.logbk$bk.tripnum.sq),a.sep="\\.",idx=1)
                 # verbose & export
                 if(any(colnames(.vms) %in% "fuelcons")){ fuelcons <- TRUE }  else fuelcons <- FALSE
                 export.check.merging.quality(a.case=1, general=general, print.console=TRUE, is.check=is.check,
                   fuelcons=fuelcons, agg.logbk.this.vessel.method.3= agg.logbk.this.vessel.method.3, .vms=.vms)
                 # find common keys
                 tripnum.sq.day.logbk            <- .logbk$bk.tripnum.sq.day
                 tripnum.sq.day.vms              <- .vms$bk.tripnum.sq.day
                 tripnum.sq.logbk                <- .logbk$bk.tripnum.sq
                 tripnum.sq.vms                  <- .vms$bk.tripnum.sq
                 tripnum.sq.day.in.vms.and.in.bk <- tripnum.sq.day.vms [tripnum.sq.day.vms %in% tripnum.sq.day.logbk]
                 tripnum.sq.in.vms.and.in.bk     <- tripnum.sq.vms [tripnum.sq.vms %in% tripnum.sq.logbk]
                 .vms.in.bk                      <- .vms[ .vms$bk.tripnum.sq.day %in%  tripnum.sq.day.in.vms.and.in.bk,]
                 .vms.in.bk2                     <- .vms[ .vms$bk.tripnum.sq %in%  tripnum.sq.in.vms.and.in.bk,]
                 if(general$landings.redistribution=="mixture123"){
                      in.bk.and.feffort.not.at.0   <- unique(.vms.in.bk[.vms.in.bk$state==1,]$bk.tripnum.sq.day)
                      in.bk2.and.feffort.not.at.0   <- unique(.vms.in.bk2[.vms.in.bk2$state==1,]$bk.tripnum.sq)
                 }
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
                                 list(.logbk.for.meth1$vesselid, .logbk.for.meth1$bk.tripnum,
                                            .logbk.for.meth1$fleet, .logbk.for.meth1$metier, .logbk.for.meth1$gear), sum, na.rm=TRUE)
                    colnames(.logbk.for.meth1) <- c("vesselid", "bk.tripnum", "fleet", "metier", "gear", nm[idx.col])
                    # verbose & export
                    export.check.merging.quality(a.case=2, general=general, print.console=TRUE,  is.check=is.check,
                                    fuelcons=fuelcons, .logbk.for.meth1=.logbk.for.meth1, .vms.for.meth1=.vms.for.meth1)
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
                                 list(.logbk.for.meth2$vesselid, .logbk.for.meth2$bk.tripnum.sq,
                                            .logbk.for.meth2$fleet, .logbk.for.meth2$metier, .logbk.for.meth2$gear), sum, na.rm=TRUE)
                    colnames(.logbk.for.meth2) <- c("vesselid", "bk.tripnum.sq", "fleet", "metier", "gear", nm[idx.col])
                    # verbose
                    export.check.merging.quality(a.case=3, general=general, print.console=TRUE,  is.check=is.check,
                       fuelcons=fuelcons, .logbk.for.meth2=.logbk.for.meth2, .vms.for.meth2=.vms.for.meth2)
                    # do.merging
                    merged2 <- do.merging(method="bk.tripnum.sq", .logbk.for.meth2, .vms.for.meth2, general)
                    # add meth flag
                    merged2$flag <- 2 # meth 2
                 }
                 #!! METH3 !!#
                 if(nrow(.logbk.for.meth3)!=0 && nrow(.vms.for.meth3)!=0 ) {
                    # verbose
                    export.check.merging.quality(a.case=4, general=general, print.console=TRUE,  is.check=is.check,
                       fuelcons=fuelcons, .logbk.for.meth3=.logbk.for.meth3, .vms.for.meth3=.vms.for.meth3)
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
                 export.check.merging.quality(a.case=5, general=general, print.console=TRUE,  is.check=is.check, merged=merged)
                 # if still 'not merging' part, retrieve on NA side i.e. occurs when pings in vms but not in bk
                   merged <- retrieve.on.bk.side(merged, type.data=c( "fleet","metier"))  # i.e. when metier=='NA'

                 # verbose
                 if(fuelcons) export.check.merging.quality(a.case=6, general=general, print.console=TRUE,  is.check=is.check,
                                     fuelcons=fuelcons,  merged=merged, inp.w=inp.w, inp.v=inp.v, inp.f=inp.f)
                 if(!fuelcons) export.check.merging.quality(a.case=6, general=general, print.console=TRUE, fuelcons=fuelcons,
                                       merged=merged, inp.w=inp.w, inp.v=inp.v)
                 export.check.merging.quality(a.case="export.this.vessel", general=general)


        # restore eflalo names
        names(merged)  [names(merged) %in% "vesselid"]   <- "VE_REF"
        names(merged)  [names(merged) %in% "fleet"]      <- "VE_FLT"
        names(merged)  [names(merged) %in% "metier"]     <- "LE_MET"
        names(merged)  [names(merged) %in% "gear"]       <- "LE_GEAR"
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
        names(merged)  [names(merged) %in% "lat"] <- "SI_LATI"
        names(merged)  [names(merged) %in% "long"] <- "SI_LONG"
        names(merged)  [names(merged) %in% "date.in.R.date"] <- "SI_DATE"
        names(merged)  [names(merged) %in% "date.in.R.time"] <- "SI_TIME"
        names(merged)  [names(merged) %in% "state"] <- "SI_STATE"
        names(merged)  [names(merged) %in% "which"] <- "SI_HARB"

       # save------------
       save("merged",   file=file.path(general$main.path,
             paste("merged.",  a.vesselid,".", general$landings.redistribution,".",general$a.year,".RData", sep='')))
       cat(paste("save 'merged'...OK\n\n",sep=""))
       } # end "mixture123"

               }else{  # end 'flag'
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



  # general settings
  general <- list()
  general$main.path    <- file.path("C:", "delivery_WP4_merging_proc_in_R")
  general$metier.def   <- "gear_targetpca" # gear_meshsize, #gear_meshsize_targetpca, or NULL if LE_MET already exists
  general$a.year       <- '2008'
  general$visual.check <- TRUE # plot for checking the first merging
  general$export.check.merging.quality <- TRUE  # export in a file the info we can see also displayed on the console...
  general$landings.redistribution <- "mixture123" # AT THE PING SCALE
   # (optional) create the export file for checking the merging quality
  general$sp.to.keep <- c("LE_KG_COD","LE_KG_PLE") 
  #=> caution: sp to keep for export check, but all species will be in 'merged' anyway.
 


  # load logbooks (eflalo)
  load(file.path(general$main.path,
     paste("eflalo.RData",sep=''))) 

  # load vms (tacsat)
  load(file.path(general$main.path,
     paste("tacsat.RData",sep=''))) 


  # TEST FOR GIVEN VESSELS
  merge.vms.to.logbook.at.the.ping.scale (logbooks=eflalo, vms=tacsat, general=general, 
                 a.vesselid=c("vessel1","vessel2"))
  #=> per vessel, merge logbook with vms
  gc(reset=TRUE)


  # read the quality check table
  qu <- read.table(file=file.path(general$main.path,
                   paste("export-check-merging-quality-eflalo-",general$a.year,".txt",sep="")), header=TRUE)
 

} # end main
