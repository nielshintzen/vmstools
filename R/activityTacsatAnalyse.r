activityTacsatAnalyse <- function(tacsat,units="year",analyse.by="LE_GEAR",identify="peaks"){

  if(!"LE_GEAR" %in% colnames(tacsat)) stop("Provide gear type (as column 'LE_GEAR' and if unknown, provide it as 'MIS'")
  if(!analyse.by %in% c("LE_GEAR","VE_REF")) warning("Analysing by unknown column variable, please check!")
  if(analyse.by %in% colnames(tacsat)){
    if(units == "all"){   yrs <- 0; mths <- 0; wks <- 0}
    if(units == "year"){  yrs <- sort(unique(format(tacsat$SI_DATIM,"%Y"))); mths  <- 0;                                    wks  <- 0}
    if(units == "month"){ yrs <- sort(unique(format(tacsat$SI_DATIM,"%Y"))); mths  <- sort(unique(month(tacsat$SI_DATIM))); wks  <- 0}
    if(units == "week"){  yrs <- sort(unique(format(tacsat$SI_DATIM,"%Y"))); wks   <- sort(unique(week(tacsat$SI_DATIM)));  mths <- 0}

    runScheme                 <- expand.grid(years=yrs,months=mths,weeks=wks,stringsAsFactors=FALSE)
    storeScheme               <- expand.grid(years=yrs,months=mths,weeks=wks,analyse.by=unique(tacsat[,analyse.by]),stringsAsFactors=FALSE)
    storeScheme$peaks         <- NA
    storeScheme$fixPeaks      <- FALSE
    storeScheme$sigma0        <- 0.911
    if(identify == "means") storeScheme$means         <- NA
    for(iRun in 1:nrow(runScheme)){
      yr  <- runScheme[iRun,"years"]
      mth <- runScheme[iRun,"months"]
      wk  <- runScheme[iRun,"weeks"]
      if(nrow(runScheme)==1){ sTacsat <- tacsat
      } else {
          if(mth == 0 & wk == 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr)
          if(mth == 0 & wk != 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & week( tacsat$SI_DATIM) == wk)
          if(mth != 0 & wk == 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & month(tacsat$SI_DATIM) == mth)
        }

      for(iBy in na.omit(unique(sTacsat[,analyse.by]))){
        dat     <- subset(sTacsat,sTacsat[,analyse.by] == iBy)
        datmr   <- dat; datmr$SI_SP <- -1*dat$SI_SP; datmr <- rbind(dat,datmr)
        xrange  <- pmin(20,range(datmr$SI_SP),na.rm=TRUE); datmr$SI_SP[which(abs(datmr$SI_SP) >20)] <- NA
        hist(datmr$SI_SP,breaks=seq(-20,20,0.5),main=paste(iBy,ifelse(yr>0,yr,""),ifelse(mth>0,mth,""),ifelse(wk>0,wk,"")),xaxt="n")
        axis(1,at=seq(-20,20,1))

        require(tcltk)

        pks <- callNumberPeak()
        storeScheme[which(storeScheme$years == yr & storeScheme$months    == mth &
                          storeScheme$weeks == wk & storeScheme$analyse.by== iBy),"peaks"] <- pks

        if(identify=="means"){
          valPeaks <- callPeakValue(pks)
          if(substr(valPeaks,1,1)==" ")
            valPeaks <- substr(valPeaks,2,nchar(valPeaks))
          storeScheme[which(storeScheme$years == yr & storeScheme$months    == mth &
                          storeScheme$weeks == wk & storeScheme$analyse.by== iBy),"means"] <- valPeaks

        }
      }
    }
  } else { stop("analyse.by statement not found as a column in the specified tacsat dataset")}
return(storeScheme)}