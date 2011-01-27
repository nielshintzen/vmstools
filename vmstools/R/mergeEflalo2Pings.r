mergeEflalo2Tacsat <-
           function(eflalo, tacsat, method="arrdeptime",...){

  lstargs <- list(...)

  #- List all vessel ID's
  all.vesselid        <- ac(unique(eflalo[anf(eflalo$VE_LEN)>=0,]$VE_REF))
  all.vesselid        <- all.vesselid[!is.na(all.vesselid)] # e.g. when VE_LEN at NA exists
  
  #- Add datim slots to both tacsat and eflalo + fix bug where time notation equals 24:00 which should be 00:00
  timeNotation        <- ifelse(length(unlist(strsplit(ac(eflalo$FT_DTIME[2]),":")))>2,"secs","mins")
  if(timeNotation == "mins"){
    res1 <- lapply(as.list(c("FT_DTIME","FT_LTIME","LE_STIME","LE_ETIME")),function(x){ idx <- grep("24:00",eflalo[,x]); eflalo[idx,x]  <- "00:00"; return(ac(eflalo[,x]))})
    res2 <- lapply(as.list(c("SI_TIME")),                                  function(x){ idx <- grep("24:00",tacsat[,x])  ; tacsat[idx,x]    <- "00:00"; return(ac(tacsat[,x]))})
    eflalo[,c("FT_DTIME","FT_LTIME","LE_STIME","LE_ETIME")]               <- do.call(cbind,res1)
    tacsat[,c("SI_TIME")]                                                 <- do.call(cbind,res2)
  } else {
      res1 <- lapply(as.list(c("FT_DTIME","FT_LTIME","LE_STIME","LE_ETIME")),function(x){ idx <- grep("24:00:00",eflalo[,x]); eflalo[idx,x]  <- "00:00"; return(eflalo[,x])})
      res2 <- lapply(as.list(c("SI_TIME")),                                  function(x){ idx <- grep("24:00:00",tacsat[,x])  ; tacsat[idx,x]    <- "00:00"; return(tacsat[,x])})
      eflalo[,c("FT_DTIME","FT_LTIME","LE_STIME","LE_ETIME")]             <- do.call(cbind,res1)
      tacsat[,c("SI_TIME")]                                               <- do.call(cbind,res2)
    }
  if(!"SI_DATIM"  %in% colnames(eflalo))  eflalo$SI_DATIM     <- as.POSIXct(paste(eflalo$FT_DDAT, eflalo$FT_DTIME,  sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
  if(!"LE_DDATIM" %in% colnames(eflalo))  eflalo$LE_DDATIM    <- eflalo$SI_DATIM
  if(!"LE_LDATIM" %in% colnames(eflalo))  eflalo$LE_LDATIM    <- as.POSIXct(paste(eflalo$FT_LDAT, eflalo$FT_LTIME,  sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
  if(!"SI_DATIM"  %in% colnames(tacsat))  tacsat$SI_DATIM     <- as.POSIXct(paste(tacsat$SI_DATE,   tacsat$SI_TIME,     sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
  
  if(length(lstargs$a.vesselid)!=0) all.vesselid <- lstargs$a.vesselid
  tacsat              <- subset(tacsat,VE_REF %in% all.vesselid)
  eflalo              <- subset(eflalo,VE_REF %in% all.vesselid)
  eflalo$VE_REF       <- factor(eflalo$VE_REF) # refactoring

  #- LOGBOOK (EFLALO) INPUT REQUIRES AT LEAST,
  #   'VE_REF',  FT_DDAT, FT_DTIME, FT_LDAT, FT_LTIME, FT_CDAT,
  #   'LE_SP_KG' (etc.), 'LE_RECT', 'VE_FLT' AND 'LE_MET_level6', 'LE_GEAR' COLUMNS
  #

  #- Tacsat: load traj with 'at sea' pings SI_STATE informed
  # ABSOLUTELY REQUIRED: c("VE_REF","SI_LATI","SI_LONG", "SI_DATE", "SI_TIME", "SI_FT", "SI_HARB", "SI_STATE")


  # keep only the essential
  if(!"SI_FT" %in% colnames(tacsat)) tacsat$SI_FT <- NA
  if(method == "midtime"){
    tacsat$SI_HARB      <- NA
    tacsat$SI_HARB      <- pointInHarbour(lon=anf(tacsat$SI_LONG), lat=anf(tacsat$SI_LATI), harbours=euharbours, rowSize=30, returnNames=F)

    # assign a trip identifier
    tacsat$SI_FT        <- 1 # init
    idx                 <- which(inHarb==0)
    tacsat[idx,"SI_FT"] <- cumsum(inHarb)[idx]
  }
  if("FT_REF" %in% colnames(tacsat)) tacsat$SI_FT <- tacsat$FT_REF
  tacsat              <- tacsat[,c("VE_REF","SI_LATI","SI_LONG",
                                                 "SI_DATIM","SI_FT")]

  # filter if vessel with a bad vms
  if(method == "midtime"){
    res                 <- by(tacsat$SI_FT,tacsat$VE_REF,FUN=function(x){length(unique(x))},simplify=T)
    tacsat              <- subset(tacsat,!VE_REF %in% names(res[which(res < 2)])) #need more than 2 trips
    idx                 <- by(tacsat$SI_FT,tacsat$VE_REF,FUN=function(x){any(is.na(x))},simplify=T)
    tacsat              <- subset(tacsat,!VE_REF %in% names(idx)[which(idx == T)]) #SI_FT is not allowed to be NA
  }

  # filter if vessel with a bad logbook
  eflalo              <- subset(eflalo,!FT_REF %in% names(which(table(eflalo$FT_REF) < 2))) #need more than 1 logbook event per vessel

  # needed if some vessels have been removed
  eflalo$VE_REF       <- factor(eflalo$VE_REF) # refactoring
  tacsat$VE_REF       <- factor(tacsat$VE_REF) # refactoring
                       
   ## ADD A WARNING IN CASE OF LONG (UNREALISTIC) TRIPS ##
  if(FALSE){
    for(i in unique(tacsat$SI_FT)){
     res <- tacsat$SI_DATIM[which(tacsat$SI_FT == i)]
     if(difftime(last(res),res[1],unit="days")>30){
        print(paste("at least one vms trip > 30 days detected! check harbours...", "\n", sep=""))
        print(paste("Check tacsat with SI_FT =",i,sep=" "))
      }
    }
  }


  if(method == "midtime"){
    # FIND THE MID-TIME OF VMS TRIPS
    startTrip           <- which(!duplicated(tacsat$SI_FT))
    endTrip             <- c((startTrip - 1)[-1],nrow(tacsat))

    res                 <- cbind(ac(unique(tacsat$SI_FT)),ac(difftime(tacsat[endTrip,"SI_DATIM"],tacsat[startTrip,"SI_DATIM"])/2 + tacsat[startTrip,"SI_DATIM"]))
    colnames(res)       <- c("SI_FT","LE_MIDTIME")
    tacsat              <- merge(tacsat,res,by="SI_FT",all=T); tacsat$LE_MIDTIME <- as.POSIXct(tacsat$LE_MIDTIME,tz="GMT",format="%d/%m/%Y  %H:%M")

    #- Calculate midpoints from eflalo: add 10 minutes and 1 minute to LDATIM and DDATIM because of bug
    eflalo$LE_MIDTIME <- difftime(eflalo$LE_LDATIM + (1*60),eflalo$LE_DDATIM+(10*60),units="secs")/2 + eflalo$LE_DDATIM + (10*60)
  }

  # THE CORE CODE: compare bk$LE_MIDTIME and vms$SI_MIDTIME
  # find the nearest bk$LE_MIDTIME for each vms$SI_MIDTIME
  # and then change levels
  # (so, for each mid.time in vms, a FT_REF will be find)
  # (so, no lines in vms without a FT_REF from bk...)

  
  eflalo              <- orderBy(~VE_REF+LE_DDATIM+LE_LDATIM,data=eflalo)
  tacsat              <- orderBy(~VE_REF+SI_DATIM,data=tacsat)
  tacefmatch          <- pmatch(sort(unique(tacsat$VE_REF)),sort(unique(eflalo$VE_REF)))

  tacsat              <- split(tacsat,tacsat$VE_REF); gc(TRUE); flush.console()
  spEflalo            <- split(eflalo,eflalo$VE_REF)

  if(method == "midtime"){
    #merge tacsat to eflalo
    for(i in 1:length(tacefmatch)){
      pm <- tacefmatch[i]
      if(is.na(pm)==T){
        tacsat[[i]]$FT_REF <- 0
      } else {
          uniqueLogMid          <- unique(spEflalo[[pm]][,c("FT_REF","LE_MIDTIME")])
          uniqueTacMid          <- unique(tacsat[[i]]$LE_MIDTIME)
          res                   <- apply(abs(outer(uniqueLogMid$LE_MIDTIME,uniqueTacMid,difftime)),2,which.min)
          tacsat[[i]]$FT_REF  <- rep(uniqueLogMid[res,"FT_REF"],table(tacsat[[i]]$LE_MIDTIME))
        }
    }
    FT_REFS                     <- unlist(lapply(tacsat,function(x){return(x$FT_REF)}))

    #non-merged eflalo to tacsat
    nonLogFT_REF                <- unique(eflalo$FT_REF)[which(!unique(eflalo$FT_REF) %in% unique(tacsat$FT_REF))]
    for(i in 1:length(nonLogFT_REF)){
      idx                       <- which(eflalo$FT_REF == nonLogFT_REF[i])
      VE_REF                    <- ac(eflalo$VE_REF[idx][1])
      Ta                        <- tacsat[[VE_REF]]
      if(is.null(dim(Ta)[1])==F){
        eflalo$FT_REF[idx]      <- rep(Ta$FT_REF[which.min(abs(difftime(eflalo$LE_MIDTIME[idx][1],Ta$LE_MIDTIME)))],length(idx))
      }
    }


  }
  if(method == "arrdeptime"){
    for(i in 1:length(tacefmatch)){
      pm    <- tacefmatch[i]
      eftim <- spEflalo[[pm]][which(duplicated(spEflalo[[pm]]$FT_REF)==F),c("LE_DDATIM","LE_LDATIM","FT_REF")]
      dtime <- eftim[,1]
      ltime <- eftim[,2]
      stime <- tacsat[[i]]$SI_DATIM
      tripn <- eftim[,3]

      if(is.na(tacefmatch[i])==T)
      {
      tacsat[[i]]$FT_REF <- 0
      }

      else {

        smdtime <- t(outer(stime,dtime,"-"))
        gtltime <- outer(ltime,stime,"-")

        #-Find first point where tacsat time is greater or equal to departure time and smaller than arrival time
        st <- apply(smdtime,1,function(x){which(x>=0)[1]})
        en <- apply(gtltime,1,function(x){rev(which(x>=0))[1]})

        #-Make sure that values are within the interval of departure and arrival time
        subse <- which(is.na(st <= en) == F & (st <= en) == T)

        st <- st[subse]
        en <- en[subse]

        #-Assign Tacsat data with FT_REF from Eflalo2 dataset where they link

        if(length(st)!=1){

          idx   <- unlist(mapply(seq,st,en,SIMPLIFY=FALSE))
          reps  <- unlist(lapply(mapply(seq,st,en,SIMPLIFY=FALSE),length))
          tacsat[[i]]$FT_REF      <- 0
          tacsat[[i]]$FT_REF[idx] <- rep(tripn[subse],reps)
        }
        if(length(st)==1){
          tacsat[[i]]$FT_REF <- 0
          tacsat[[i]]$FT_REF[seq(st,en)] <- rep(tripn[subse],length(seq(st,en)))
        }
        if(length(st)==0){

          tacsat[[i]]$FT_REF <- 0
        }
      }
    }

    FT_REFS     <- unlist(lapply(tacsat,function(x){return(x$FT_REF)}))
  }
  
  tacsat        <- do.call(rbind,tacsat)
  tacsat$FT_REF <- FT_REFS

 
return(list(eflalo=eflalo,tacsat=tacsat))}