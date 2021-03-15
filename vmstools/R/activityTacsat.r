#' Define activity of fishing vessel based on speed profile analyses.
#' 
#' Analyse tacsat data by gear or vessel and define, based on speed profile and
#' fitting normal distributions through these speed profiles, fishing and other
#' activities.
#' 
#' The analysis assumes that the speed profile close to 0 can be mirrored to
#' create a good normal distribution. Success of fit might depend on the
#' initial guess in peaks or mean peak values.
#' 
#' Note that the default value of sigma0 = 0.911 (used when prior estimates of
#' parameters are not given). The value of 0.911 corresponds to a width of
#' 1.5knots on each side of the mean under a 90 percent CI.
#' 
#' @param tacsat tacsat dataset
#' @param units Analyse by: "year", "month" and "week". "month" and "week"
#' cannot be used at same time.
#' @param analyse.by Analyse tacsat by gear ("LE_GEAR") or vessel ("VE_REF").
#' @param storeScheme If \code{\link{activityTacsatAnalyse}} is used, supply
#' output here.
#' @param plot Logical. Whether the results of the fit of the normal
#' distributions should be plotted
#' @param level If analyse.by="LE_GEAR" is selected, there is an option to
#' analyse at vessel level too, but taking the gear parameter estimates as
#' starting values. level = "vessel" turns this option on.
#' @return In general, a 5-peak analysis results in h=no fishing / in harbour,
#' f=fishing, s=steaming. A 3-peak analyses results in f=fishing (closest to
#' zero speed) and s=steaming.
#' @author Niels T. Hintzen
#' @seealso \code{\link{activityTacsatAnalyse}},
#' \code{\link{segmentedTacsatSpeed}}
#' @examples
#' 
#' data(tacsat)
#' data(eflalo)
#' 
#' tacsatp         <- mergeEflalo2Tacsat(eflalo,tacsat)
#' tacsatp$LE_GEAR <- eflalo$LE_GEAR[match(tacsatp$FT_REF,eflalo$FT_REF)]
#' 
#' tacsatp$LE_GEAR <- ac(tacsatp$LE_GEAR)
#' tacsatp$LE_GEAR[which(is.na(tacsatp$LE_GEAR)==TRUE)] <- "NO_GEAR"
#' tacsat          <- tacsatp
#' 
#' #--------- LE_GEAR -----------
#' #- Visual analyses of activity, and define number of peaks / kernels that can be
#' #  observed (either 3 or 5)
#' \dontrun{
#' storeScheme   <- activityTacsatAnalyse(subset(tacsat,LE_GEAR == "OTM" &
#'                       format(SI_DATIM,"%Y") == 1801),units="year",analyse.by="LE_GEAR",identify="means")
#' res           <- activityTacsat(subset(tacsat,LE_GEAR == "OTM" &
#'                       format(SI_DATIM,"%Y") == 1801),units="year",analyse.by="LE_GEAR",
#'                       storeScheme,plot=TRUE)
#' res           <- activityTacsat(subset(tacsat,LE_GEAR == "OTM" &
#'                       format(SI_DATIM,"%Y") == 1801),units="year",analyse.by="LE_GEAR",
#'                       storeScheme,plot=TRUE,level="vessel")
#' 
#' #--------- VE_REF -----------
#' tacsat        <- subset(tacsat,VE_REF == "801")
#' storeScheme   <- activityTacsatAnalyse(tacsat,units="year",analyse.by="VE_REF",
#'                                        identify="means")
#' 
#' #- Run activityTacsat
#' res           <- activityTacsat(tacsat,units="year",analyse.by="VE_REF",storeScheme,
#'                                 plot=TRUE,level="all")
#' }
#' 
#' @export activityTacsat
activityTacsat <- function(tacsat,units="year",analyse.by="LE_GEAR",storeScheme=NULL,plot=FALSE,level="all"){

  require("mixtools")
  if (!"SI_DATIM" %in% colnames(tacsat))
        tacsat$SI_DATIM <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,
            sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M")
  if(!analyse.by %in% c("LE_GEAR","VE_REF")) stop("Analysing only by gear or vessel")

  #- Make subset for only those tacsat records that have speed
  tacsat$ID   <- 1:nrow(tacsat)
  tacsat$SI_STATE <- NA
  tacsatOrig  <- tacsat
  idx         <- which(is.na(tacsat$SI_SP)==FALSE)
  tacsat      <- tacsat[idx,]
  
  #- If sigma is NULL it needs to be estimated and gets a variable name
  storeScheme$sigma0[which(storeScheme$sigma0==0)] <- "d"

  if(units == "all"){   yrs <- 0;                                   mths  <- 0;                                    wks  <- 0}
  if(units == "year"){  yrs <- sort(unique(format(tacsat$SI_DATIM,"%Y"))); mths  <- 0;                                    wks  <- 0}
  if(units == "month"){ yrs <- sort(unique(format(tacsat$SI_DATIM,"%Y"))); mths  <- an(sort(unique(format(tacsat$SI_DATIM,"%m")))); wks  <- 0}
  if(units == "week"){  yrs <- sort(unique(format(tacsat$SI_DATIM,"%Y"))); wks   <- an(sort(unique(format(tacsat$SI_DATIM,"%W"))))+1;  mths <- 0}
  
  runScheme                 <- expand.grid(years=yrs,months=mths,weeks=wks)
  
  #-----------------------------------------------------------------------------
  # Start run for all combinations of gear / vessel and units
  #-----------------------------------------------------------------------------
  
  for(iRun in 1:nrow(runScheme)){
    yr  <- runScheme[iRun,"years"]
    mth <- runScheme[iRun,"months"]
    wk  <- runScheme[iRun,"weeks"]
    if(nrow(runScheme)==1){ sTacsat <- tacsat
    } else {
        if(mth == 0 & wk == 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr)
        if(mth == 0 & wk != 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & (an(format( tacsat$SI_DATIM,"%W"))+1) == wk)
        if(mth != 0 & wk == 0) sTacsat <- subset(tacsat,format(tacsat$SI_DATIM,"%Y") == yr & format(tacsat$SI_DATIM,"%m") == mth)
      }
    if(plot==TRUE) x11();
    #---------------------------------------------------------------------------
    #- Analyses when gear info is supplied  LE_GEAR
    #---------------------------------------------------------------------------
    if("LE_GEAR" %in% colnames(tacsat) & analyse.by == "LE_GEAR"){

      gearList  <- names(which((rowSums(table(sTacsat$LE_GEAR,sTacsat$SI_SP)) - table(sTacsat$LE_GEAR,sTacsat$SI_SP)[,"0"])>40))

      #- Mirror the tacsat dataset and make a selection
      tyg       <- subset(sTacsat,is.na(LE_GEAR) == FALSE  & LE_GEAR %in% gearList); tygmr <- tyg; tygmr$SI_SP <- -1* tygmr$SI_SP; tygmr <- rbind(tyg,tygmr)
      tng       <- subset(sTacsat,is.na(LE_GEAR) == TRUE | !LE_GEAR %in% gearList); tngmr <- tng; tngmr$SI_SP <- -1* tngmr$SI_SP; tngmr <- rbind(tng,tngmr)

      #-------------------------------------------------------------------------
      #- Get general speed pattern by gear, use analysed number of kernals LE_GEAR + GENERIC
      #-------------------------------------------------------------------------
      res       <- list()
      for(iGr in unique(tyg$LE_GEAR)){

        #- Get rid of very influential datapoints (lower their abundance)
          tbl                 <- table(subset(tygmr,LE_GEAR==iGr)$SI_SP);
          spd                 <- an(names(rev(sort(tbl))[1]))
          idx                 <- which(subset(tygmr,LE_GEAR==iGr)$SI_SP==spd)
          nxt                 <- ifelse(names(rev(sort(tbl))[1])==ac(spd),ifelse(abs(an(names(rev(sort(tbl))[2])))==abs(spd),names(rev(sort(tbl))[3]),names(rev(sort(tbl))[2])),names(rev(sort(tbl))[1]))
          if(tbl[ac(spd)]/tbl[nxt] > 5){
            idx <- sample(idx,tbl[ac(spd)]-tbl[nxt]*2,replace=FALSE)
            if(length(which(abs(an(names(tbl))) %in% spd))>1) idx <- c(idx,sample(which(subset(tygmr,LE_GEAR==iGr)$SI_SP==(-1*spd)),tbl[ac(-1*spd)]-tbl[nxt]*2,replace=FALSE))
          } else { idx <- -1:-nrow(subset(tygmr,LE_GEAR==iGr))}

        #-----------------------------------------------------------------------------
        # Fit the 3 or 5 normal distributions. If parameter guestimates are
        #  available, then use these
        #-----------------------------------------------------------------------------
        if(is.null(storeScheme)==TRUE){
          res[[iGr]]    <- try(normalmixEM(subset(tygmr,LE_GEAR==iGr)$SI_SP[-idx],maxit=1000,k=5,maxrestarts=20,mean.constr=c("-b","-a",0,"a","b"),sd.constr=c("b","a",0.911,"a","b"),sigma=rep(1,5)))
        } else {
            #- Fitting model when mean values of peaks has been defined
            if("means" %in% colnames(storeScheme)){

              #- Extract parameters from storeScheme
              ss          <- storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                               storeScheme$weeks == wk & storeScheme$analyse.by == iGr),"means"]
              sigma       <- anf(storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                               storeScheme$weeks == wk & storeScheme$analyse.by == iGr),"sigma0"])
              fixPeaks    <- ac( storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                               storeScheme$weeks == wk & storeScheme$analyse.by == iGr),"fixPeaks"])

              #- Setup parameter estimate vectors for mu and sigma
              if(length(c(na.omit(as.numeric(strsplit(ss," ")[[1]]))))==3){ constraintmn <- c("-a",0,"a") } else { constraintmn <- c("-b","-a",0,"a","b")}
              if(length(c(na.omit(as.numeric(strsplit(ss," ")[[1]]))))==3){ constraintsd <- c("a","b","a")} else { constraintsd <- c("b","a",sigma,"a","b")}
              if(fixPeaks) constraintmn <- c(na.omit(anf(unlist(strsplit(ss," ")))))

              #- Fit the actual model through the normalmixEM function
              res[[iGr]]   <- try(normalmixEM(subset(tygmr,LE_GEAR==iGr)$SI_SP[-idx],maxit=1000,mu=c(na.omit(as.numeric(strsplit(ss," ")[[1]]))), sigma=rep(1,length(constraintsd)),
                                              maxrestarts=20,mean.constr=constraintmn,sd.constr=constraintsd))
            } else {
              #- Fitting model when number of peaks has been defined
              
              #- Extract parameters from storeScheme
              ss          <- storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                               storeScheme$weeks == wk & storeScheme$analyse.by == iGr),"peaks"]
              sigma       <- anf(storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                               storeScheme$weeks == wk & storeScheme$analyse.by == iGr),"sigma0"])
              fixPeaks    <- ac( storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                               storeScheme$weeks == wk & storeScheme$analyse.by == iGr),"fixPeaks"])

              #- Setup parameter estimate vectors for mu and sigma
              if(ss==3){ constraintmn <- c("-a",0,"a") } else { constraintmn <- c("-b","-a",0,"a","b")}
              if(ss==3){ constraintsd <- c("a","b","a")} else { constraintsd <- c("b","a",sigma,"a","b")}
              if(length(ss)>0){

                #- Fit the actual model through the normalmixEM function
                if(is.na(ss)==TRUE)  res[[iGr]]    <- try(normalmixEM(subset(tygmr,LE_GEAR==iGr)$SI_SP[-idx],maxit=1000,k=5, maxrestarts=20,  mean.constr=c("-b","-a",0,"a","b"),sd.constr=c("b","a",sigma,"a","b"),sigma=rep(1,5)))
                if(is.na(ss)==FALSE) res[[iGr]]    <- try(normalmixEM(subset(tygmr,LE_GEAR==iGr)$SI_SP[-idx],maxit=1000,k=ss,maxrestarts=20,  mean.constr=constraintmn,          sd.constr=constraintsd,            sigma=rep(1,length(constraintsd))))
              } else {               res[[iGr]]    <- try(normalmixEM(subset(tygmr,LE_GEAR==iGr)$SI_SP[-idx],maxit=1000,k=5, maxrestarts=20,  mean.constr=c("-b","-a",0,"a","b"),sd.constr=c("b","a",sigma,"a","b"),sigma=rep(1,5)))}
            }
          }
        if(plot==TRUE) plot(res[[iGr]],2,breaks=100,xlim=c(-20,20))
      }
      if(level == "vessel"){
        #- Transform the output into the right format
        for(iGr in unique(tyg$LE_GEAR))
          if(!class(res[[iGr]]) == "try-error"){ res[[iGr]] <- res[[iGr]]$mu } else { res[[iGr]] <- rep(NA,5)}
        res             <- lapply(res,function(x){if(class(x)=="try-error"){x<-rep(NA,5)}else{x}})
        res             <- lapply(res,sort)
      }

      if(level == "vessel"){
        #-------------------------------------------------------------------------
        #- Perform analyses per vessel with gear LE_GEAR + VE_REF
        #-------------------------------------------------------------------------

        if(nrow(tygmr)>40)
          shipList  <- names(which((rowSums(table(tygmr$VE_REF,tygmr$SI_SP)) - table(tygmr$VE_REF,tygmr$SI_SP)[,"0"])>20))
        shipFit   <- list()
        if(exists("shipList")){
          for(iShip in shipList){

            #- Get rid of very influential datapoints (lower their abundance)
            tbl                 <- table(subset(tygmr,VE_REF==iShip)$SI_SP);
            spd                 <- an(names(rev(sort(tbl))[1]))
            idx                 <- which(subset(tygmr,VE_REF==iShip)$SI_SP==spd)
            nxt                 <- ifelse(names(rev(sort(tbl))[1])==ac(spd),ifelse(abs(an(names(rev(sort(tbl))[2])))==abs(spd),names(rev(sort(tbl))[3]),names(rev(sort(tbl))[2])),names(rev(sort(tbl))[1]))
            if(tbl[ac(spd)]/tbl[nxt] > 5){
              idx <- sample(idx,tbl[ac(spd)]-tbl[nxt]*2,replace=FALSE)
              if(length(which(abs(an(names(tbl))) %in% spd))>1) idx <- c(idx,sample(which(subset(tygmr,VE_REF==iShip)$SI_SP==(-1*spd)),tbl[ac(-1*spd)]-tbl[nxt]*2,replace=FALSE))
            } else { idx <- -1:-nrow(subset(tygmr,VE_REF==iShip))}

            shipTacsat        <- subset(tygmr,VE_REF == iShip)

            #-----------------------------------------------------------------------------
            # Fit the 3 or 5 normal distributions. If parameter guestimates are
            #  available, then use these
            #-----------------------------------------------------------------------------

            #- Setup parameter estimate vectors for mu and sigma
            if(length(res[[names(which.max(table(shipTacsat$LE_GEAR)))]])==3){ constraintmn <- c("-a",0,"a")} else { constraintmn <- c("-b","-a",0,"a","b")}
            if(length(res[[names(which.max(table(shipTacsat$LE_GEAR)))]])==3){ constraintsd <- c("a","b","a")}else { constraintsd <- c("b","a",0.911,"a","b")}

            #- Fit the actual model through the normalmixEM function
            shipFit[[iShip]]  <- try(normalmixEM(shipTacsat$SI_SP[-idx],mu=res[[names(which.max(table(shipTacsat$LE_GEAR)))]],maxit=2000,
                                                 sigma=rep(1,length(constraintsd)),mean.constr=constraintmn,sd.constr=constraintsd))

            if(class(shipFit[[iShip]])!= "try-error"){

              #- Analyse the fit and turn it into a result of fishing - no fishing
              #mu                <- sort.int(shipFit[[iShip]]$mu,index.return=TRUE)
              #sds               <- shipFit[[iShip]]$sigma[mu$ix]; mu <- mu$x
              mu                <- shipFit[[iShip]]$mu
              sds               <- shipFit[[iShip]]$sigma

              probs             <- dnorm(x=shipTacsat$SI_SP,mean=mu[ceiling(length(mu)/2)],sd=sds[ceiling(length(mu)/2)])
              for(i in (ceiling(length(mu)/2)+1):length(mu)) probs <- cbind(probs,dnorm(x=shipTacsat$SI_SP,mean=mu[i],sd=sds[i]))
              SI_STATE          <- apply(probs,1,which.max)
              
              if(length(mu)==3){
                SI_STATE        <- af(SI_STATE); levels(SI_STATE) <- c("f","s"); SI_STATE <- ac(SI_STATE)}
              if(length(mu)==5){
                SI_STATE        <- af(SI_STATE); levels(SI_STATE) <- c("h","f","s"); SI_STATE <- ac(SI_STATE)}
              tacsat$SI_STATE[which(tacsat$ID %in% shipTacsat$ID)] <- SI_STATE[1:(length(SI_STATE)/2)]
            } else { tacsat$SI_STATE[which(tacsat$ID %in% shipTacsat$ID)] <- NA}
          }
        }
        } else {
          for(iGr in unique(tyg$LE_GEAR)){
            if(!class(res[[iGr]]) == "try-error"){

              #- Analyse the fit and turn it into a result of fishing - no fishing
              #mu                  <- sort.int(res[[iGr]]$mu,index.return=TRUE)
              #sds                 <- res[[iGr]]$sigma[mu$ix]; mu <- mu$x
              mu                  <- res[[iGr]]$mu
              sds                 <- res[[iGr]]$sigma
              probs               <- dnorm(x=subset(tyg,LE_GEAR==iGr)$SI_SP,mean=mu[ceiling(length(mu)/2)],sd=sds[ceiling(length(mu)/2)])
              for(i in (ceiling(length(mu)/2)+1):length(mu)) probs <- cbind(probs,dnorm(x=subset(tyg,LE_GEAR==iGr)$SI_SP,mean=mu[i],sd=sds[i]))
              SI_STATE            <- apply(probs,1,which.max)

              if(length(mu)==3){
                SI_STATE        <- af(SI_STATE); levels(SI_STATE) <- c("f","s"); SI_STATE <- ac(SI_STATE)}
              if(length(mu)==5){
                SI_STATE        <- af(SI_STATE); levels(SI_STATE) <- c("h","f","s"); SI_STATE <- ac(SI_STATE)}
              tacsat$SI_STATE[which(tacsat$ID %in% subset(tyg,LE_GEAR == iGr)$ID)] <- SI_STATE
            }
          }
        }
        #-------------------------------------------------------------------------
        #- Perform analyses per vessel without gear NO_GEAR + VE_REF
        #-------------------------------------------------------------------------
        if(nrow(tngmr)>40)
          nonshipList           <- names(which((rowSums(table(tngmr$VE_REF,tngmr$SI_SP)) - table(tngmr$VE_REF,tngmr$SI_SP)[,"0"])>20))
        nonshipFit            <- list()
        if(exists("nonshipList")){
          for(iShip in nonshipList){

            #- Get rid of very influential datapoints (lower their abundance)
            tbl                 <- table(subset(tngmr,VE_REF==iShip)$SI_SP);
            spd                 <- an(names(rev(sort(tbl))[1]))
            idx                 <- which(subset(tngmr,VE_REF==iShip)$SI_SP==spd)
            nxt                 <- ifelse(names(rev(sort(tbl))[1])==ac(spd),ifelse(abs(an(names(rev(sort(tbl))[2])))==abs(spd),names(rev(sort(tbl))[3]),names(rev(sort(tbl))[2])),names(rev(sort(tbl))[1]))
            if(tbl[ac(spd)]/tbl[nxt] > 5){
              idx <- sample(idx,tbl[ac(spd)]-tbl[nxt]*2,replace=FALSE)
              if(length(which(abs(an(names(tbl))) %in% spd))>1) idx <- c(idx,sample(which(subset(tngmr,VE_REF==iShip)$SI_SP==(-1*spd)),tbl[ac(-1*spd)]-tbl[nxt]*2,replace=FALSE))
            } else { idx <- -1:-nrow(subset(tngmr,VE_REF==iShip))}

            #-----------------------------------------------------------------------------
            # Fit the 3 or 5 normal distributions. If parameter guestimates are
            #  available, then use these
            #-----------------------------------------------------------------------------

            shipTacsat          <- subset(tngmr,VE_REF == iShip)
            if(exists("shipFit")){
              if(iShip %in% names(shipFit)){

                #- Setup parameter estimate vectors for mu and sigma
                if(length(shipFit[[iShip]]$mu)==3){constraintmn <- c("-a",0,"a")} else { constraintmn <- c("-b","-a",0,"a","b")}
                if(length(shipFit[[iShip]]$mu)==3){constraintsd <- c("a","b","a")}else { constraintsd <- c("b","a",0.911,"a","b")}

                #- Fit the actual model through the normalmixEM function
                nonshipFit[[iShip]] <- try(normalmixEM(shipTacsat$SI_SP[-idx],k=length(shipFit[[iShip]]$mu),maxit=2000,
                                                       sigma=rep(1,length(constraintsd)),mean.constr=constraintmn,sd.constr=constraintsd))
            } else {
                #- Fit the actual model through the normalmixEM function
                nonshipFit[[iShip]] <- try(normalmixEM(shipTacsat$SI_SP[-idx],k=5,maxit=2000,mean.constr=c("-b","-a",0,"a","b"),sd.constr=c("b","a",0.911,"a","b"),sigma=rep(1,5)))}

            if(!class(nonshipFit[[iShip]]) == "try-error"){

              #- Analyse the fit and turn it into a result of fishing - no fishing
              #mu                  <- sort.int(nonshipFit[[iShip]]$mu,index.return=TRUE)
              #sds                 <- nonshipFit[[iShip]]$sigma[mu$ix]; mu <- mu$x
              mu                  <- nonshipFit[[iShip]]$mu
              sds                 <- nonshipFit[[iShip]]$sds

              probs               <- dnorm(x=shipTacsat$SI_SP,mean=mu[ceiling(length(mu)/2)],sd=sds[ceiling(length(mu)/2)])
              for(i in (ceiling(length(mu)/2)+1):length(mu)) probs <- cbind(probs,dnorm(x=shipTacsat$SI_SP,mean=mu[i],sd=sds[i]))
              SI_STATE            <- apply(probs,1,which.max)
              
              if(length(mu)==3){
                SI_STATE        <- af(SI_STATE); levels(SI_STATE) <- c("f","s"); SI_STATE <- ac(SI_STATE)}
              if(length(mu)==5){
                SI_STATE        <- af(SI_STATE); levels(SI_STATE) <- c("h","f","s"); SI_STATE <- ac(SI_STATE)}
              tacsat$SI_STATE[which(tacsat$ID %in% shipTacsat$ID)] <- SI_STATE[1:(length(SI_STATE)/2)]
            }
          }
        }
      }
    }
    #---------------------------------------------------------------------------
    #- Analyses when vessel info is supplied  VE_REF + ENOUGH #
    #---------------------------------------------------------------------------
    if("VE_REF" %in% colnames(tacsat) & analyse.by == "VE_REF"){

      vesselList  <- names(which((rowSums(table(sTacsat$VE_REF,sTacsat$SI_SP)) - table(sTacsat$VE_REF,sTacsat$SI_SP)[,"0"])>40))

      #- Mirror the tacsat dataset and make a selection
      tyv       <- subset(sTacsat,is.na(VE_REF) == FALSE  & VE_REF %in% vesselList); tyvmr <- tyv; tyvmr$SI_SP <- -1* tyvmr$SI_SP; tyvmr <- rbind(tyv,tyvmr)
      tnv       <- subset(sTacsat,is.na(VE_REF) == TRUE | !VE_REF %in% vesselList); tnvmr <- tnv; tnvmr$SI_SP <- -1* tnvmr$SI_SP; tnvmr <- rbind(tnv,tnvmr)

      #- Perform analyses per vessel
      if(nrow(tyv)>40)
        shipList  <- names(which((rowSums(table(tyvmr$VE_REF,tyvmr$SI_SP)) - table(tyvmr$VE_REF,tyvmr$SI_SP)[,"0"])>20))
      shipFit   <- list()
      if(exists("shipList")){
        for(iShip in shipList){

          #- Get rid of very influential data points
          tbl                 <- table(subset(tyvmr,VE_REF==iShip)$SI_SP);
          spd                 <- an(names(rev(sort(tbl))[1]))
          idx                 <- which(subset(tyvmr,VE_REF==iShip)$SI_SP==spd)
          nxt                 <- ifelse(names(rev(sort(tbl))[1])==ac(spd),ifelse(abs(an(names(rev(sort(tbl))[2])))==abs(spd),names(rev(sort(tbl))[3]),names(rev(sort(tbl))[2])),names(rev(sort(tbl))[1]))
          if(tbl[ac(spd)]/tbl[nxt] > 5){
            idx <- sample(idx,tbl[ac(spd)]-tbl[nxt]*2,replace=FALSE)
            if(length(which(abs(an(names(tbl))) %in% spd))>1) idx <- c(idx,sample(which(subset(tyvmr,VE_REF==iShip)$SI_SP==(-1*spd)),tbl[ac(-1*spd)]-tbl[nxt]*2,replace=FALSE))
          } else { idx <- -1:-nrow(subset(tyvmr,VE_REF==iShip))}

          shipTacsat        <- subset(tyvmr,VE_REF == iShip)

          #-----------------------------------------------------------------------------
          # Fit the 3 or 5 normal distributions. If parameter guestimates are
          #  available, then use these
          #-----------------------------------------------------------------------------
          if(is.null(storeScheme)==TRUE){
            shipFit[[iShip]]    <- try(normalmixEM(subset(tyvmr,VE_REF==iShip)$SI_SP[-idx],maxit=2000,k=5,mean.constr=c("-b","-a",0,"a","b"),sd.constr=c("b","a",0.911,"a","b"),sigma=rep(1,5)))
          } else {
              #- Fitting model when mean values of peaks has been defined
              if("means" %in% colnames(storeScheme)){

                #- Extract parameters from storeScheme
                ss          <- storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                                 storeScheme$weeks == wk & storeScheme$analyse.by == iShip),"means"]
                sigma       <- anf(storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                                 storeScheme$weeks == wk & storeScheme$analyse.by == iShip),"sigma0"])
                fixPeaks    <- ac( storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                                 storeScheme$weeks == wk & storeScheme$analyse.by == iShip),"fixPeaks"])

                #- Setup parameter estimate vectors for mu and sigma
                if(length(c(na.omit(as.numeric(strsplit(ss," ")[[1]]))))==3){ constraintmn <- c("-a",0,"a")} else { constraintmn <- c("-b","-a",0,"a","b")}
                if(length(c(na.omit(as.numeric(strsplit(ss," ")[[1]]))))==3){ constraintsd <- c("a","b","a")}else { constraintsd <- c("b","a",sigma,"a","b")}
                if(fixPeaks) constraintmn <- c(na.omit(anf(unlist(strsplit(ss," ")))))

                #- Fit the actual model through the normalmixEM function
                shipFit[[iShip]]   <- try(normalmixEM(subset(tyvmr,VE_REF==iShip)$SI_SP[-idx],mu=c(na.omit(as.numeric(strsplit(ss," ")[[1]]))), maxit=2000,
                                           mean.constr=constraintmn,sd.constr=constraintsd,sigma=rep(1,length(constraintsd))))
              } else {
                #- Fitting model when number of peaks has been defined

                #- Extract parameters from storeScheme
                ss          <- storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                                 storeScheme$weeks == wk & storeScheme$analyse.by == iShip),"peaks"]
                sigma       <- anf(storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                                 storeScheme$weeks == wk & storeScheme$analyse.by == iShip),"sigma0"])
                fixPeaks    <- ac( storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                                 storeScheme$weeks == wk & storeScheme$analyse.by == iShip),"fixPeaks"])

                #- Setup parameter estimate vectors for mu and sigma
                if(ss==3){ constraintmn <- c("-a",0,"a")} else { constraintmn <- c("-b","-a",0,"a","b")}
                if(ss==3){ constraintsd <- c("a","b","a")}else { constraintsd <- c("b","a",sigma,"a","b")}
                if(length(ss)>0){

                  #- Fit the actual model through the normalmixEM function
                  if(is.na(ss)==TRUE) shipFit[[iShip]]   <- try(normalmixEM(subset(tyvmr,VE_REF==iShip)$SI_SP[-idx],maxit=2000,k=5,mean.constr=c("-b","-a",0,"a","b"),sd.constr=c("b","a",sigma,"a","b"),sigma=rep(1,5)))
                  if(is.na(ss)==FALSE) shipFit[[iShip]]   <- try(normalmixEM(subset(tyvmr,VE_REF==iShip)$SI_SP[-idx],maxit=2000,k=ss,mean.constr=constraintmn,sd.constr=constraintsd,sigma=rep(1,length(constraintsd))))
                } else {           shipFit[[iShip]]   <- try(normalmixEM(subset(tyvmr,VE_REF==iShip)$SI_SP[-idx],maxit=2000,k=5,mean.constr=c("-b","-a",0,"a","b"),sd.constr=c("b","a",sigma,"a","b"),sigma=rep(1,5)))}
              }
            }
          if(plot==TRUE) plot(shipFit[[iShip]],2,breaks=100,xlim=c(-20,20))
          if(!class(shipFit[[iShip]]) == "try-error"){

            #- Analyse the fit and turn it into a result of fishing - no fishing
            #mu                <- sort.int(shipFit[[iShip]]$mu,index.return=TRUE)
            #sds               <- shipFit[[iShip]]$sigma[mu$ix]; mu <- mu$x
            mu                <- shipFit[[iShip]]$mu
            sds               <- shipFit[[iShip]]$sigma

            probs             <- dnorm(x=shipTacsat$SI_SP,mean=mu[ceiling(length(mu)/2)],sd=sds[ceiling(length(mu)/2)])
            for(i in (ceiling(length(mu)/2)+1):length(mu)) probs <- cbind(probs,dnorm(x=shipTacsat$SI_SP,mean=mu[i],sd=sds[i]))
            SI_STATE          <- apply(probs,1,which.max)

            if(length(mu)==3){
              SI_STATE        <- af(SI_STATE); levels(SI_STATE) <- c("f","s"); SI_STATE <- ac(SI_STATE)}
            if(length(mu)==5){
              SI_STATE        <- af(SI_STATE); levels(SI_STATE) <- c("h","f","s"); SI_STATE <- ac(SI_STATE)}
            tacsat$SI_STATE[which(tacsat$ID %in% shipTacsat$ID)] <- SI_STATE[1:(length(SI_STATE)/2)]
          }
        }
      }
      #-------------------------------------------------------------------------
      #- Perform analyses for all vessels left over VE_REF + NOT ENOUGH #
      #-------------------------------------------------------------------------
      if(nrow(tnvmr)>40)
        nonshipList           <- names(which((rowSums(table(tnvmr$VE_REF,tnvmr$SI_SP)) - table(tnvmr$VE_REF,tnvmr$SI_SP)[,"0"])>20))
      nonshipFit            <- list()
      if(exists("nonshipList")){
        for(iShip in nonshipList){

          #- Get rid of very influential datapoints (lower their abundance)
          tbl                 <- table(subset(tnvmr,VE_REF==iShip)$SI_SP);
          spd                 <- an(names(rev(sort(tbl))[1]))
          idx                 <- which(subset(tnvmr,VE_REF==iShip)$SI_SP==spd)
          nxt                 <- ifelse(names(rev(sort(tbl))[1])==ac(spd),ifelse(abs(an(names(rev(sort(tbl))[2])))==abs(spd),names(rev(sort(tbl))[3]),names(rev(sort(tbl))[2])),names(rev(sort(tbl))[1]))
          if(tbl[ac(spd)]/tbl[nxt] > 5){
            idx <- sample(idx,tbl[ac(spd)]-tbl[nxt]*2,replace=FALSE)
            if(length(which(abs(an(names(tbl))) %in% spd))>1) idx <- c(idx,sample(which(subset(tnvmr,VE_REF==iShip)$SI_SP==(-1*spd)),tbl[ac(-1*spd)]-tbl[nxt]*2,replace=FALSE))
          } else { idx <- -1:-nrow(subset(tnvmr,VE_REF==iShip))}

          shipTacsat          <- subset(tnvmr,VE_REF == iShip)
          #-----------------------------------------------------------------------------
          # Fit the 3 or 5 normal distributions. If parameter guestimates are
          #  available, then use these
          #-----------------------------------------------------------------------------
          if(length(shipFit[[iShip]]$mu)==3){constraintmn <- c("-a",0,"a")} else { constraintmn <- c("-b","-a",0,"a","b")}
          if(length(shipFit[[iShip]]$mu)==3){constraintsd <- c("a","b","a")}else { constraintsd <- c("b","a",0.911,"a","b")}

          #- Fit the actual model through the normalmixEM function
          nonshipFit[[iShip]] <- try(normalmixEM(shipTacsat$SI_SP[-idx],k=length(shipFit[[iShip]]$mu),maxit=2000,mean.constr=constraintmn,sd.constr=constraintsd,sigma=rep(1,length(constraintsd))))

          if(!class(nonshipFit[[iShip]]) == "try-error"){

            #- Analyse the fit and turn it into a result of fishing - no fishing
            #mu                  <- sort.int(nonshipFit[[iShip]]$mu,index.return=TRUE)
            #sds                 <- nonshipFit[[iShip]]$sigma[mu$ix]; mu <- mu$x
            mu                  <- nonshipFit[[iShip]]$mu
            sds                 <- nonshipFit[[iShip]]$sigma

            probs               <- dnorm(x=shipTacsat$SI_SP,mean=mu[ceiling(length(mu)/2)],sd=sds[ceiling(length(mu)/2)])
            for(i in (ceiling(length(mu)/2)+1):length(mu)) probs <- cbind(probs,dnorm(x=shipTacsat$SI_SP,mean=mu[i],sd=sds[i]))
            SI_STATE            <- apply(probs,1,which.max)

            if(length(mu)==3)
              SI_STATE        <- af(SI_STATE); levels(SI_STATE) <- c("f","s"); SI_STATE <- ac(SI_STATE)
            if(length(mu)==5)
              SI_STATE        <- af(SI_STATE); levels(SI_STATE) <- c("h","f","s"); SI_STATE <- ac(SI_STATE)
            tacsat$SI_STATE[which(tacsat$ID %in% shipTacsat$ID)] <- SI_STATE[1:(length(SI_STATE)/2)]
          }
        }
      }
    }
  }
  
leftOverTacsat  <- tacsatOrig[which(!tacsatOrig$ID %in% tacsat$ID),]
tacsat          <- rbind(tacsat,leftOverTacsat)
tacsat          <- orderBy(~ID,tacsat)
cat("Note that in case of 5 peaks: no fishing = h, fishing = f, steaming / no fishing = s\n")
cat("Note that in case of 3 peaks: fishing = f, steaming / no fishing = s\n")
return(tacsat$SI_STATE)}
