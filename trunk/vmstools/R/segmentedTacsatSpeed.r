segmentedTacsatSpeed <- function(tacsat,units="year",analyse.by="VE_REF",speed="calculated",logfit=FALSE,CI=0.95,saveDir=tempdir(),forceLowerBound=NULL){

  require(segmented)
  tacsat$idxFun <- 1:nrow(tacsat)

  #Regular checks
  if(!'SI_STATE' %in% colnames(tacsat)) tacsat$SI_STATE  <- NA
  if(!"SI_DATIM" %in% colnames(tacsat)) tacsat$SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")
  if(analyse.by == "LE_GEAR"){ if(!"LE_GEAR" %in% colnames(tacsat)) stop("Provide gear type (as column 'LE_GEAR' and if unknown, provide it as 'MIS'")}
  if(!analyse.by %in% c("LE_GEAR","VE_REF","VE_REF+LE_GEAR")) warning("Analysing by unknown column variable, please check!")
  if(analyse.by == "VE_REF+LE_GEAR"){
    if(!"VE_REF" %in% colnames(tacsat) | !"LE_GEAR" %in% colnames(tacsat)) stop("VE_REF and LE_GEAR not both in tacsat available")
    tacsat$VE_REF_LE_GEAR <- paste(tacsat$VE_REF,tacsat$LE_GEAR)
    analyse.by <- "VE_REF_LE_GEAR"
  }

  tacsatOrig <- tacsat
  if(analyse.by %in% colnames(tacsat)){
    if(units == "all"){   yrs <- 0; mths <- 0; wks <- 0}
    if(units == "year"){  yrs <- sort(unique(format(tacsat$SI_DATIM,"%Y"))); mths  <- 0;                                    wks  <- 0}
    if(units == "month"){ yrs <- sort(unique(format(tacsat$SI_DATIM,"%Y"))); mths  <- sort(unique(month(tacsat$SI_DATIM))); wks  <- 0}
    if(units == "week"){  yrs <- sort(unique(format(tacsat$SI_DATIM,"%Y"))); wks   <- sort(unique(week(tacsat$SI_DATIM)));  mths <- 0}
  } else { stop("analyse.by statement not found as a column in the specified tacsat dataset")}
    storeScheme               <- expand.grid(years=yrs,months=mths,weeks=wks,analyse.by=unique(tacsat[,analyse.by]))
    # Add upper and lower boundaries to storeScheme
    storeScheme               <- cbind(storeScheme,data.frame(lower=rep(0,nrow(storeScheme)),
                                                              upper=rep(0,nrow(storeScheme))))
    storeScheme$analyse.by    <- ac(storeScheme$analyse.by)
    storeScheme$success       <- 0
  for(iRun in 1:nrow(storeScheme)){
      yr  <- storeScheme[iRun,"years"]
      mth <- storeScheme[iRun,"months"]
      wk  <- storeScheme[iRun,"weeks"]
      aby <- storeScheme[iRun,"analyse.by"]
      if(analyse.by == "VE_REF"){
        if(nrow(storeScheme)==1){ sTacsat <- tacsat
        } else {
            if(mth == 0 & wk == 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & VE_REF == aby)
            if(mth == 0 & wk != 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & week( tacsat$SI_DATIM) == wk & VE_REF == aby)
            if(mth != 0 & wk == 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & month(tacsat$SI_DATIM) == mth & VE_REF == aby)
          }
      }
      if(analyse.by == "LE_GEAR"){
        if(nrow(storeScheme)==1){ sTacsat <- tacsat
        } else {
            if(mth == 0 & wk == 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & LE_GEAR == aby)
            if(mth == 0 & wk != 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & week( tacsat$SI_DATIM) == wk & LE_GEAR == aby)
            if(mth != 0 & wk == 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & month(tacsat$SI_DATIM) == mth & LE_GEAR == aby)
          }
      }
      if(analyse.by == "VE_REF_LE_GEAR"){
        if(nrow(storeScheme)==1){ sTacsat <- tacsat
        } else {
            if(mth == 0 & wk == 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & VE_REF_LE_GEAR == aby )
            if(mth == 0 & wk != 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & week( tacsat$SI_DATIM) == wk & VE_REF_LE_GEAR == aby)
            if(mth != 0 & wk == 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & month(tacsat$SI_DATIM) == mth & VE_REF_LE_GEAR == aby)
          }
      }
      
      #Check if there is any data available
      if(nrow(sTacsat)>0){

        #Check if you want to run with calculated or instantaneous speed
        if(speed == "calculated"){
          sTacsat             <- sortTacsat(sTacsat)
          sTacsat$SI_SP_ORIG  <- sTacsat$SI_SP
          if("FT_REF" %in% colnames(tacsat)){ sTacsat$SI_SP <- calculateSpeed(sTacsat,level="trip",   weight=c(0.5,0.5), fill.na=TRUE)$SI_SPCA
          } else {                             sTacsat$SI_SP <- calculateSpeed(sTacsat,level="vessel", weight=c(0.5,0.5), fill.na=TRUE)$SI_SPCA}
        }
        #Remove records where SI_SP is NA
        sTacsat <- sTacsat[which(is.na(sTacsat$SI_SP)==F & sTacsat$SI_SP > 0 & sTacsat$SI_SP <= 10),]


        minRows <- 20 #Based on a scan of sub tacsat datasets, 5% percentile was approx 20

        #Exception exercise: If number of rows in sTacsat is really low, return NA
        if(nrow(sTacsat) <= minRows){
          storeScheme[iRun,"lower"] <- NA
          storeScheme[iRun,"upper"] <- NA
        }

        #Regular case: If number of rows is enough to fit segmented regression
        if(nrow(sTacsat) > minRows){

          #Define starting values for segmented regression
          hi  <- hist(sTacsat$SI_SP,breaks=diff(c(floor(  range(sTacsat$SI_SP[is.finite(sTacsat$SI_SP)],na.rm=TRUE)[1]),
                                                  ceiling(range(sTacsat$SI_SP[is.finite(sTacsat$SI_SP)],na.rm=TRUE)[2]))),plot=FALSE)
          acc <- diff(diff(cumsum(hi$counts)))
          idx <- rev(sort(abs(diff(diff(cumsum(hi$counts))))))[1:2]
          cnts<- which(abs(acc) %in% idx)+2 #Taking twice diff, so add 2 to get back to counts
          psi <- list(x=range(hi$breaks[cnts]))               #First guess on breakpoints
          psiOrig <- psi
          dat <- data.frame(x=sort(sTacsat$SI_SP[sTacsat$SI_SP>0]),
                            y=1:length(sTacsat$SI_SP[sTacsat$SI_SP>0]))
          if(logfit==TRUE) dat <- data.frame(x=      an(rep(names(table(sTacsat$SI_SP[sTacsat$SI_SP>0])),
                                                         ceiling(log(table(sTacsat$SI_SP[sTacsat$SI_SP>0]))))),
                                          y=1:length(rep(names(table(sTacsat$SI_SP[sTacsat$SI_SP>0])),
                                                         ceiling(log(table(sTacsat$SI_SP[sTacsat$SI_SP>0]))))))
          o <- 1 ; class(o) <- "try-error" ; count <- 0  ; bound1 <-NULL ; bound2 <- NULL;

          #Fit the model
          while(class(o)=="try-error"){
            count <- count+1
            o <- try(
                     segmented(lm(y~x, data=dat) , seg.Z=~x , psi=psi, control= seg.control(display = FALSE, it.max=50, h=1)), # with 2 starting guesses
                     silent=TRUE) # the second breakpoint is quite uncertain and could lead to failure so...
            if(!"try-error" %in% class(o)) break else psi <- list(x=c(sort(runif(2,min=range(sTacsat$SI_SP,na.rm=TRUE)[1],max=range(sTacsat$SI_SP,na.rm=TRUE)[2])))) # searching decreasing by 1 each time
            if(count>20) {bound1 <- psiOrig$x[1]; bound2 <- psiOrig$x[2] ; cat("failure of the segmented regression for",paste(c("year","month","week","analyse.by"),storeScheme[iRun,1:4]),"\n"); break}
          }
          #Calculate the bounds and whether the fit has been successful or not
          if(is.null(bound1)==T & is.null(bound2)==TRUE){
            bound1 <- max(range(sTacsat$SI_SP)[1],min(confint(o,level=CI)$x))
            bound2 <- min(range(sTacsat$SI_SP)[2],max(confint(o,level=CI)$x))
            if(class(o)[1] != "try-error") storeScheme[iRun,"success"] <- 1
          }

          #Save the bounds
          if(is.null(forceLowerBound)==FALSE)
            bound1 <- forceLowerBound
          if(bound2 < bound1)
            bound2 <- bound1
          storeScheme[iRun,"lower"] <- bound1
          storeScheme[iRun,"upper"] <- bound2

          tacsatOrig$SI_STATE[sTacsat$idxFun[which(sTacsat$SI_SP >= storeScheme[iRun,"lower"] & sTacsat$SI_SP <= storeScheme[iRun,"upper"])]] <- "f" #Fishing
          tacsatOrig$SI_STATE[sTacsat$idxFun[which(sTacsat$SI_SP <  storeScheme[iRun,"lower"] | sTacsat$SI_SP >  storeScheme[iRun,"upper"])]] <- "nf" #Steaming / in harbour
        }
      }
    }
  #Write the results to file and display the success rates
  write.csv(storeScheme,file=file.path(saveDir,"storeScheme.csv"))
  cat("Successful segmented regression fits",length(which(storeScheme$success==1)),"\n",
      "versus unsuccessful fits",length(which(storeScheme$success == 0)),"\n\n",
      "Check ",file.path(saveDir,"storeScheme.csv"),"for details \n\n")
      
  cat("Note: fishing = f, no fishing = nf\n")
return(tacsatOrig[,-grep("idxFun",colnames(tacsatOrig))])}
    
#res <- segmentedTacsatSpeed(tacsat,units="year",analyse.by="VE_REF",speed="calculated",logfit=FALSE,CI=0.95)
    

