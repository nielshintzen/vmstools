activityTacsat <- function(tacsat,units="year",analyse.by="LE_GEAR",storeScheme=NULL){

  require("mixtools")
  if (!"SI_DATIM" %in% colnames(tacsat))
        tacsat$SI_DATIM <- as.POSIXct(paste(tacsat$SI_DATE, tacsat$SI_TIME,
            sep = " "), tz = "GMT", format = "%d/%m/%Y  %H:%M")
  if(!analyse.by %in% c("LE_GEAR","VE_REF")) stop("Analysing only by gear or vessel")

  #- Make subset for only those tacsat records that have speed
  tacsat$ID   <- 1:nrow(tacsat)
  tacsat$SI_STATE <- NA
  tacsatOrig  <- tacsat
  idx         <- which(is.na(tacsat$SI_SP)==F)
  tacsat      <- tacsat[idx,]

  if(units == "all"){   yrs <- 0;                                   mths  <- 0;                                    wks  <- 0}
  if(units == "year"){  yrs <- sort(unique(year(tacsat$SI_DATIM))); mths  <- 0;                                    wks  <- 0}
  if(units == "month"){ yrs <- sort(unique(year(tacsat$SI_DATIM))); mths  <- sort(unique(month(tacsat$SI_DATIM))); wks  <- 0}
  if(units == "week"){  yrs <- sort(unique(year(tacsat$SI_DATIM))); wks   <- sort(unique(week(tacsat$SI_DATIM)));  mths <- 0}
  
  runScheme                 <- expand.grid(years=yrs,months=mths,weeks=wks)
  
  #- BUG fix to normalmix.init as embedded within mixtools
  normalmix.init <- function (x, lambda = NULL, mu = NULL, s = NULL, k = 2, arbmean = TRUE,arbvar = TRUE){
    if (!is.null(s)) {arbvar <- (length(s) > 1)
        if (arbvar) k <- length(s)
    }
    if (!is.null(mu)) {arbmean <- (length(mu) > 1)
        if (arbmean) { k <- length(mu)
            if (!is.null(s) && length(s) > 1 && k != length(s)) {
                stop("mu and sigma are each of length >1 but not of the same length.")}}}
    if (!arbmean && !arbvar) { stop("arbmean and arbvar cannot both be FALSE")}
    n = length(x);x = sort(x);x.bin = list()
    for (j in 1:k) {
      x.bin[[j]] <- x[max(1, floor((j - 1) * n/k)):ceiling(j * n/k)]}
    if (is.null(s)) {
        s.hyp = as.vector(sapply(x.bin, sd))
        s.hyp[which(s.hyp==0)] <- runif(sum(s.hyp==0),0,sd(x))
        if (arbvar) { s = 1/rexp(k, rate = s.hyp)
        } else { s = 1/rexp(1, rate = mean(s.hyp))}
    }
    if (is.null(mu)) {
        mu.hyp <- as.vector(sapply(x.bin, mean))
        if (arbmean) { mu = rnorm(k, mean = mu.hyp, sd = s)
        } else { mu = rnorm(1, mean = mean(mu.hyp), sd = mean(s))}
    }
    if (is.null(lambda)) {lambda <- runif(k);lambda <- lambda/sum(lambda)
    } else {lambda <- rep(lambda, length.out = k);lambda <- lambda/sum(lambda)}
    list(lambda = lambda, mu = mu, s = s, k = k, arbvar = arbvar,
        arbmean = arbmean)
  }
  
  #-----------------------------------------------------------------------------
  # Start run for all combinations of gear / vessel and units
  #-----------------------------------------------------------------------------
  
  for(iRun in 1:nrow(runScheme)){
    yr  <- runScheme[iRun,"years"]
    mth <- runScheme[iRun,"months"]
    wk  <- runScheme[iRun,"weeks"]
    if(nrow(runScheme)==1){ sTacsat <- tacsat
    } else {
        if(mth == 0 & wk == 0) sTacsat <- subset(tacsat,year(tacsat$SI_DATIM) == yr)
        if(mth == 0 & wk != 0) sTacsat <- subset(tacsat,year(tacsat$SI_DATIM) == yr & week( tacsat$SI_DATIM) == wk)
        if(mth != 0 & wk == 0) sTacsat <- subset(tacsat,year(tacsat$SI_DATIM) == yr & month(tacsat$SI_DATIM) == mth)
      }
    #---------------------------------------------------------------------------
    #- Analyses when gear info is supplied  LE_GEAR
    #---------------------------------------------------------------------------
    if("LE_GEAR" %in% colnames(tacsat) & analyse.by == "LE_GEAR"){

      gearList  <- names(which((rowSums(table(sTacsat$LE_GEAR,sTacsat$SI_SP)) - table(sTacsat$LE_GEAR,sTacsat$SI_SP)[,"0"])>40))

      #- Mirror the tacsat dataset and make a selection
      tyg       <- subset(sTacsat,is.na(LE_GEAR) == F  & LE_GEAR %in% gearList); tygmr <- tyg; tygmr$SI_SP <- -1* tygmr$SI_SP; tygmr <- rbind(tyg,tygmr)
      tng       <- subset(sTacsat,is.na(LE_GEAR) == T | !LE_GEAR %in% gearList); tngmr <- tng; tngmr$SI_SP <- -1* tngmr$SI_SP; tngmr <- rbind(tng,tngmr)

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
            idx <- sample(idx,tbl[ac(spd)]-tbl[nxt]*2,replace=F)
            if(length(which(abs(an(names(tbl))) %in% spd))>1) idx <- c(idx,sample(which(subset(tygmr,LE_GEAR==iGr)$SI_SP==(-1*spd)),tbl[ac(-1*spd)]-tbl[nxt]*2,replace=F))
          } else { idx <- -1:-nrow(subset(tygmr,LE_GEAR==iGr))}

        if(is.null(storeScheme)==T){
          res[[iGr]]    <- try(normalmixEM(subset(tygmr,LE_GEAR==iGr)$SI_SP[-idx],maxit=1000,k=5,maxrestarts=10)$mu)
        } else {
            ss          <- storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                             storeScheme$weeks == wk & storeScheme$analyse.by == iGr),"peaks"]
            if(length(ss)>0){
              if(is.na(ss)==T)  res[[iGr]]   <- try(normalmixEM(subset(tygmr,LE_GEAR==iGr)$SI_SP[-idx],maxit=1000,k=5, maxrestarts=10)$mu)
              if(is.na(ss)==F)  res[[iGr]]   <- try(normalmixEM(x=subset(tygmr,LE_GEAR==iGr)$SI_SP[-idx],maxit=1000,k=ss,maxrestarts=10)$mu)
            } else {            res[[iGr]]   <- try(normalmixEM(subset(tygmr,LE_GEAR==iGr)$SI_SP[-idx],maxit=1000,k=5, maxrestarts=10)$mu)}
          }
      }
      res             <- lapply(res,function(x){if(class(x)=="try-error"){x<-rep(NA,5)}else{x}})
      res             <- lapply(res,sort)

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
            idx <- sample(idx,tbl[ac(spd)]-tbl[nxt]*2,replace=F)
            if(length(which(abs(an(names(tbl))) %in% spd))>1) idx <- c(idx,sample(which(subset(tygmr,VE_REF==iShip)$SI_SP==(-1*spd)),tbl[ac(-1*spd)]-tbl[nxt]*2,replace=F))
          } else { idx <- -1:-nrow(subset(tygmr,VE_REF==iShip))}

          shipTacsat        <- subset(tygmr,VE_REF == iShip)
          shipFit[[iShip]]  <- try(normalmixEM(shipTacsat$SI_SP[-idx],mu=res[[names(which.max(table(shipTacsat$LE_GEAR)))]],maxit=2000))

          if(class(shipFit[[iShip]])!= "try-error"){
            mu                <- sort.int(shipFit[[iShip]]$mu,index.return=T)
            sds               <- shipFit[[iShip]]$sigma[mu$ix]; mu <- mu$x

            probs             <- dnorm(x=shipTacsat$SI_SP,mean=mu[ceiling(length(mu)/2)],sd=sds[ceiling(length(mu)/2)])
            for(i in (ceiling(length(mu)/2)+1):length(mu)) probs <- cbind(probs,dnorm(x=shipTacsat$SI_SP,mean=mu[i],sd=sds[i]))
            SI_STATE          <- apply(probs,1,which.max)
            tacsat$SI_STATE[which(tacsat$ID %in% shipTacsat$ID)] <- SI_STATE
          } else { tacsat$SI_STATE[which(tacsat$ID %in% shipTacsat$ID)] <- NA}
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
            idx <- sample(idx,tbl[ac(spd)]-tbl[nxt]*2,replace=F)
            if(length(which(abs(an(names(tbl))) %in% spd))>1) idx <- c(idx,sample(which(subset(tngmr,VE_REF==iShip)$SI_SP==(-1*spd)),tbl[ac(-1*spd)]-tbl[nxt]*2,replace=F))
          } else { idx <- -1:-nrow(subset(tngmr,VE_REF==iShip))}

          shipTacsat          <- subset(tngmr,VE_REF == iShip)
          if(iShip %in% names(shipFit)) nonshipFit[[iShip]] <- try(normalmixEM(shipTacsat$SI_SP[-idx],k=length(shipFit[[iShip]]$mu),maxit=2000))
          if(!iShip%in% names(shipFit)) nonshipFit[[iShip]] <- try(normalmixEM(shipTacsat$SI_SP[-idx],k=5,maxit=2000))

          mu                  <- sort.int(nonshipFit[[iShip]]$mu,index.return=T)
          sds                 <- nonshipFit[[iShip]]$sigma[mu$ix]; mu <- mu$x

          probs               <- dnorm(x=shipTacsat$SI_SP,mean=mu[ceiling(length(mu)/2)],sd=sds[ceiling(length(mu)/2)])
          for(i in (ceiling(length(mu)/2)+1):length(mu)) probs <- cbind(probs,dnorm(x=shipTacsat$SI_SP,mean=mu[i],sd=sds[i]))
          SI_STATE            <- apply(probs,1,which.max)
          tacsat$SI_STATE[which(tacsat$ID %in% shipTacsat$ID)] <- SI_STATE
        }
      }
    }
    #---------------------------------------------------------------------------
    #- Analyses when vessel info is supplied  VE_REF + ENOUGH #
    #---------------------------------------------------------------------------
    if("VE_REF" %in% colnames(tacsat) & analyse.by == "VE_REF"){

      vesselList  <- names(which((rowSums(table(sTacsat$VE_REF,sTacsat$SI_SP)) - table(sTacsat$VE_REF,sTacsat$SI_SP)[,"0"])>40))

      #- Mirror the tacsat dataset and make a selection
      tyv       <- subset(sTacsat,is.na(VE_REF) == F  & VE_REF %in% vesselList); tyvmr <- tyv; tyvmr$SI_SP <- -1* tyvmr$SI_SP; tyvmr <- rbind(tyv,tyvmr)
      tnv       <- subset(sTacsat,is.na(VE_REF) == T | !VE_REF %in% vesselList); tnvmr <- tnv; tnvmr$SI_SP <- -1* tnvmr$SI_SP; tnvmr <- rbind(tnv,tnvmr)

      #- Perform analyses per vessel
      if(nrow(tyv)>40)
        shipList  <- names(which((rowSums(table(tyvmr$VE_REF,tyvmr$SI_SP)) - table(tyvmr$VE_REF,tyvmr$SI_SP)[,"0"])>20))
      shipFit   <- list()
      if(exists("shipList")){
        for(iShip in shipList){

          tbl                 <- table(subset(tyvmr,VE_REF==iShip)$SI_SP);
          spd                 <- an(names(rev(sort(tbl))[1]))
          idx                 <- which(subset(tyvmr,VE_REF==iShip)$SI_SP==spd)
          nxt                 <- ifelse(names(rev(sort(tbl))[1])==ac(spd),ifelse(abs(an(names(rev(sort(tbl))[2])))==abs(spd),names(rev(sort(tbl))[3]),names(rev(sort(tbl))[2])),names(rev(sort(tbl))[1]))
          if(tbl[ac(spd)]/tbl[nxt] > 5){
            idx <- sample(idx,tbl[ac(spd)]-tbl[nxt]*2,replace=F)
            if(length(which(abs(an(names(tbl))) %in% spd))>1) idx <- c(idx,sample(which(subset(tyvmr,VE_REF==iShip)$SI_SP==(-1*spd)),tbl[ac(-1*spd)]-tbl[nxt]*2,replace=F))
          } else { idx <- -1:-nrow(subset(tyvmr,VE_REF==iShip))}

          shipTacsat        <- subset(tyvmr,VE_REF == iShip)

          if(is.null(storeScheme)==T){
            shipFit[[iShip]]    <- try(normalmixEM(subset(tyvmr,VE_REF==iShip)$SI_SP[-idx],maxit=2000,k=5)$mu)
          } else {
              ss          <- storeScheme[which(storeScheme$years == yr & storeScheme$months     == mth &
                                               storeScheme$weeks == wk & storeScheme$analyse.by == iShip),"peaks"]
              if(length(ss)>0){
                if(is.na(ss)==T) shipFit[[iShip]]   <- try(normalmixEM(subset(tyvmr,VE_REF==iShip)$SI_SP[-idx],maxit=2000,k=5))
                if(is.na(ss)==F) shipFit[[iShip]]   <- try(normalmixEM(subset(tyvmr,VE_REF==iShip)$SI_SP[-idx],maxit=2000,k=ss))
              } else {           shipFit[[iShip]]   <- try(normalmixEM(subset(tyvmr,VE_REF==iShip)$SI_SP[-idx],maxit=2000,k=5))}
            }

          mu                <- sort.int(shipFit[[iShip]]$mu,index.return=T)
          sds               <- shipFit[[iShip]]$sigma[mu$ix]; mu <- mu$x

          probs             <- dnorm(x=shipTacsat$SI_SP,mean=mu[ceiling(length(mu)/2)],sd=sds[ceiling(length(mu)/2)])
          for(i in (ceiling(length(mu)/2)+1):length(mu)) probs <- cbind(probs,dnorm(x=shipTacsat$SI_SP,mean=mu[i],sd=sds[i]))
          SI_STATE          <- apply(probs,1,which.max)
          tacsat$SI_STATE[which(tacsat$ID %in% shipTacsat$ID)] <- SI_STATE
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
            idx <- sample(idx,tbl[ac(spd)]-tbl[nxt]*2,replace=F)
            if(length(which(abs(an(names(tbl))) %in% spd))>1) idx <- c(idx,sample(which(subset(tnvmr,VE_REF==iShip)$SI_SP==(-1*spd)),tbl[ac(-1*spd)]-tbl[nxt]*2,replace=F))
          } else { idx <- -1:-nrow(subset(tnvmr,VE_REF==iShip))}

          shipTacsat          <- subset(tnvmr,VE_REF == iShip)
          nonshipFit[[iShip]] <- try(normalmixEM(shipTacsat$SI_SP[-idx],k=length(shipFit[[iShip]]$mu),maxit=2000))

          mu                  <- sort.int(nonshipFit[[iShip]]$mu,index.return=T)
          sds                 <- nonshipFit[[iShip]]$sigma[mu$ix]; mu <- mu$x

          probs               <- dnorm(x=shipTacsat$SI_SP,mean=mu[ceiling(length(mu)/2)],sd=sds[ceiling(length(mu)/2)])
          for(i in (ceiling(length(mu)/2)+1):length(mu)) probs <- cbind(probs,dnorm(x=shipTacsat$SI_SP,mean=mu[i],sd=sds[i]))
          SI_STATE            <- apply(probs,1,which.max)
          tacsat$SI_STATE[which(tacsat$ID %in% shipTacsat$ID)] <- SI_STATE
        }
      }
    }
  }
return(tacsat$SI_STATE)}