# do the distinction between fishing and non-fishing
# by automatic detection of the fishing peak
# adapted from library 'segmented'
# F. Bastardie

segmentTacsatSpeed <- function(tacsat,
                               vessels= unique(tacsat$VE_REF),
                               force.lower.bound=0.5,  
                               gears.to.force= c('GNS','GND','GNC','GNF','GTR','GTN','GEN','GN','SDN','SSC'),
                               general=list(
                                   output.path=file.path('C:','output'),
                                   visual.check=TRUE, 
                                   a.year=2009,
                                   what.speed="calculated", 
                                    ),
                                ...){

 
   lstargs <-  as.list( sys.call() ) # deprecated

  # checks
  if(!'SI_STATE' %in% colnames(tacsat)) tacsat$SI_STATE <- NA
  if(!'idx' %in% colnames(tacsat))  tacsat$idx <- 1:nrow(tacsat)
  cat("if it still doesn't exist, 'results' folder is created in ",general$output.path,"\n")
  dir.create(general$output.path, showWarnings = TRUE, recursive = TRUE, mode = "0777")
 
  # add
  tacsat$bound1 <- NA
  tacsat$bound2 <- NA


   # utils---
  distAB.f <- function(A,B, .unit="km"){
  # return the dist in km or nautical miles.
  # coord input: e.g. A[54,10] means [54 N, 10 E]
  # formula = cos (gamma) = sin LatA sin LatB + cos LatA cos LatB cos delta
  # if gamma in degree, knowing that 1 degree = 60 minutes = 60 nm, then gamma in degree * 60 = dist in nm
  # we can also use the radius R of the earth: the length of the arc AB is given by R*gamma with gamma in radians
        p = 180/pi;  Rearth = 6378.388 # in km
        nm <- 1852e-3 # 1852 m for 1 nm
        res <- switch(.unit,
           km = Rearth * acos(sin(A[,1]/p)*sin(B[,1]/p) + (cos(A[,1]/p) * cos(B[,1]/p)*cos(A[,2]/p - B[,2]/p))),
           nm = Rearth * acos(sin(A[,1]/p)*sin(B[,1]/p) + (cos(A[,1]/p) * cos(B[,1]/p)*cos(A[,2]/p - B[,2]/p))) / nm
           )
       names(res) <- NULL
        return(res)
        }

    # utils---
    # library(segmented)
    # from library(segmented) # getS3method("segmented", "lm")
    ## CAUTION: UNFORTUNATELY, THE FUNCTION NEED TO BE PLACED HERE BECAUSE A PROBLEM OF ENV
    seg.control <- function (toll = 1e-04, it.max = 20, display = FALSE, last = TRUE,
    maxit.glm = 25, h = 1){
    list(toll = toll, it.max = it.max, visual = display, last = last,
        maxit.glm = maxit.glm, h = h)
    }
    segmented <- function (obj, seg.Z, psi, control = seg.control(), model.frame = TRUE,
    ...){
    it.max <- old.it.max <- control$it.max
    toll <- control$toll
    visual <- control$visual
    last <- control$last
    h <- min(abs(control$h), 1)
    if (h < 1)
        it.max <- it.max + round(it.max/2)
    objframe <- update(obj, model = TRUE, x = TRUE, y = TRUE)
    y <- objframe$y
    a <- model.matrix(seg.Z, data = eval(obj$call$data))
    a <- subset(a, select = colnames(a)[-1])
    n <- nrow(a)
    Z <- list()
    for (i in colnames(a)) Z[[length(Z) + 1]] <- a[, i]
    name.Z <- names(Z) <- colnames(a)
    if (length(Z) == 1 && is.vector(psi) && is.numeric(psi)) {
        psi <- list(as.numeric(psi))
        names(psi) <- name.Z
    }
    if (!is.list(Z) || !is.list(psi) || is.null(names(Z)) ||
        is.null(names(psi)))
        stop("Z and psi have to be *named* list")
    nomiZpsi <- match(names(Z), names(psi))
    if (!identical(length(Z), length(psi)) || any(is.na(nomiZpsi)))
        stop("Length or names of Z and psi do not match")
    dd <- match(names(Z), names(psi))
    nome <- names(psi)[dd]
    psi <- psi[nome]
    a <- sapply(psi, length)
    b <- rep(1:length(a), times = a)
    Znew <- list()
    for (i in 1:length(psi)) Znew[[length(Znew) + 1]] <- rep(Z[i],
        a[i])
    Z <- matrix(unlist(Znew), nrow = n)
    colnames(Z) <- rep(nome, a)
    psi <- unlist(psi)
    k <- ncol(Z)
    PSI <- matrix(rep(psi, rep(n, k)), ncol = k)
    nomiZ <- rep(nome, times = a)
    ripetizioni <- as.numeric(unlist(sapply(table(nomiZ)[order(unique(nomiZ))],
        function(xxx) {
            1:xxx
        })))
    nomiU <- paste("U", ripetizioni, sep = "")
    nomiU <- paste(nomiU, nomiZ, sep = ".")
    KK <- new.env()
    for (i in 1:ncol(objframe$model)) assign(names(objframe$model[i]),
        objframe$model[[i]], envir = KK)
    if (it.max == 0) {
        U <- pmax((Z - PSI), 0)
        colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
        nomiU <- paste("U", colnames(U), sep = "")
        for (i in 1:ncol(U)) assign(nomiU[i], U[, i], envir = KK)
        Fo <- update.formula(formula(obj), as.formula(paste(".~.+",
            paste(nomiU, collapse = "+"))))
        obj <- update(obj, formula = Fo, data = KK)
        if (model.frame)
            obj$model <- data.frame(as.list(KK))
        obj$psi <- psi
        return(obj)
    }
    XREG <- model.matrix(obj)
    o <- model.offset(objframe)
    w <- model.weights(objframe)
    if (is.null(w))
        w <- rep(1, n)
    if (is.null(o))
        o <- rep(0, n)
    initial <- psi
    it <- 1
    epsilon <- 10
    obj0 <- obj
    list.obj <- list(obj)
    psi.values <- NULL
    rangeZ <- apply(Z, 2, range)
    H <- 1
    while (abs(epsilon) > toll) {
        U <- pmax((Z - PSI), 0)
        V <- ifelse((Z > PSI), -1, 0)
        dev.old <- sum(obj$residuals^2)

        X <- cbind(XREG, U, V)
        rownames(X) <- NULL
        if (ncol(V) == 1)
            colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c("U", "V")
        else colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c(paste("U",
            1:k, sep = ""), paste("V", 1:k, sep = ""))
        obj <- lm.wfit(x = X, y = y, w = w, offset = o)
        dev.new <- sum(obj$residuals^2)
        if (visual) {
            if (it == 1)
                cat(0, " ", formatC(dev.old, 3, format = "f"),
                  "", "(No breakpoint(s))", "\n")
            spp <- if (it < 10)
                ""
            else NULL
            cat(it, spp, "", formatC(dev.new, 3, format = "f"),
                "\n")
        }
        epsilon <- (dev.new - dev.old)/dev.old
        obj$epsilon <- epsilon
        it <- it + 1
        obj$it <- it
        class(obj) <- c("segmented", class(obj))
        list.obj[[length(list.obj) + ifelse(last == TRUE, 0,
            1)]] <- obj
        if (k == 1) {
            beta.c <- coef(obj)["U"]
            gamma.c <- coef(obj)["V"]
        }
        else {
            beta.c <- coef(obj)[paste("U", 1:k, sep = "")]
            gamma.c <- coef(obj)[paste("V", 1:k, sep = "")]
        }
        if (it > it.max)
            break
        psi.values[[length(psi.values) + 1]] <- psi.old <- psi
        if (it >= old.it.max && h < 1)
            H <- h
        psi <- psi.old + H * gamma.c/beta.c
        PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
        a <- apply((Z < PSI), 2, all)
        b <- apply((Z > PSI), 2, all)
        if (sum(a + b) != 0 || is.na(sum(a + b)))
            stop("(Some) estimated psi out of its range")
        obj$psi <- psi
    }
    psi.values[[length(psi.values) + 1]] <- psi
    id.warn <- FALSE
    if (it > it.max) {
        warning("max number of iterations attained", call. = FALSE)
        id.warn <- TRUE
    }
    Vxb <- V %*% diag(beta.c, ncol = length(beta.c))
    colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
    colnames(Vxb) <- paste(ripetizioni, nomiZ, sep = ".")
    nomiU <- paste("U", colnames(U), sep = "")
    nomiVxb <- paste("psi", colnames(Vxb), sep = "")
    for (i in 1:ncol(U)) {
        assign(nomiU[i], U[, i], envir = KK)
        assign(nomiVxb[i], Vxb[, i], envir = KK)
    }
    nnomi <- c(nomiU, nomiVxb)
    Fo <- update.formula(formula(obj0), as.formula(paste(".~.+",
        paste(nnomi, collapse = "+"))))
    objF <- update(obj0, formula = Fo, data = KK)
    Cov <- vcov(objF)
    id <- match(paste("psi", colnames(Vxb), sep = ""), names(coef(objF)))
    vv <- if (length(id) == 1)
        Cov[id, id]
    else diag(Cov[id, id])
    psi <- cbind(initial, psi, sqrt(vv))
    rownames(psi) <- colnames(Cov)[id]
    colnames(psi) <- c("Initial", "Est.", "St.Err")
    objF$rangeZ <- rangeZ
    objF$psi.history <- psi.values
    objF$psi <- psi
    objF$it <- (it - 1)
    objF$epsilon <- epsilon
    objF$call <- match.call()
    objF$nameUV <- list(U = paste("U", colnames(Vxb), sep = ""),
        V = rownames(psi), Z = name.Z)
    objF$id.warn <- id.warn
    if (model.frame)
        objF$mframe <- data.frame(as.list(KK))
    class(objF) <- c("segmented", class(obj0))
    list.obj[[length(list.obj) + 1]] <- objF
    class(list.obj) <- "segmented"
    if (last)
        list.obj <- list.obj[[length(list.obj)]]
    return(list.obj)
    }


  # hereafter, the core code...
  # A FOR-LOOP BY VESSEL--------------------
  for(a.vesselid in vessels){ 
    tacsat.this.vessel <- tacsat[tacsat$VE_REF %in% a.vesselid, ]
  
    if('SI_DATE' %in% colnames(tacsat.this.vessel))
          {
          tacsat.this.vessel$SI_DATIM  <- strptime(  paste(tacsat.this.vessel$SI_DATE, tacsat.this.vessel$SI_TIME),
                            tz='GMT', "%e/%m/%Y %H:%M" )
          } else{
          if(!any(colnames(tacsat.this.vessel) %in% 'SI_DATIM')) stop('you need either to inform a SI_DATIM or a SI_DATE')
          }
            
    if(is.null(general$speed)) general$speed <- "calculated"
    if(general$speed=="calculated"){
     # add a apparent speed colunm (nautical miles per hour)
     diff.time <- tacsat.this.vessel[-nrow(tacsat.this.vessel),"SI_DATIM"] -
                       tacsat.this.vessel[-1,"SI_DATIM"]
     tacsat.this.vessel$diff.time.mins <- c(0, as.numeric(diff.time, units="mins"))
     tacsat.this.vessel$speed <- signif(
            c(0, distAB.f(A= tacsat.this.vessel[-nrow(tacsat.this.vessel),c("SI_LATI","SI_LONG")],
              B= tacsat.this.vessel[-1,c("SI_LATI","SI_LONG")], .unit="nm") / (abs(tacsat.this.vessel[-1,]$diff.time.mins)/60)),
              3)
     }         
     if(general$speed=="observed"){
       tacsat.this.vessel$speed <- tacsat.this.vessel$SI_SP # instantaneous speed
       }          

    # cleaning irrealistic points
    tacsat.this.vessel$speed <-
       replace(tacsat.this.vessel$speed, is.na(tacsat.this.vessel$speed), 0)


    idx <- tacsat.this.vessel[tacsat.this.vessel$speed >= 30 |
                is.infinite(tacsat.this.vessel$speed),"idx"]
    tacsat <- tacsat[!tacsat$idx %in% idx,] # just remove!
    tacsat.this.vessel <- tacsat.this.vessel[tacsat.this.vessel$speed < 30,]
 
    # check
    if(is.null(levels(factor(tacsat.this.vessel$LE_GEAR))))
      stop('you need first to assign a gear LE_GEAR to each ping')
  

    # A FOR-LOOP BY GEAR--------------------
      for (gr in levels(factor(tacsat.this.vessel$LE_GEAR))){ 

        xxx <- tacsat.this.vessel[tacsat.this.vessel$LE_GEAR==gr,] # input

        x <- as.numeric(as.character(sort(xxx$speed))) *100   # multiply by factor 100 because integer needed
        hi <- hist(x, nclass=30,plot=FALSE) 

        y <- c(1:length(sort(xxx$speed))) # sort in increasing order
        y <- y[x>100 & x<1000] # just removing the 0, and the larger speeds we 100% know it is steaming
        x <- x[x>100 & x<1000]
        dati   <- data.frame(x=x,y=y)
        dati$x <- as.integer(dati$x) # integer needed
        psi    <- list(x= quantile(dati$x,probs=c(0.05,0.5))  )
        assign('dati', dati, envir=.GlobalEnv) # DEBUG segmented()...this function looks in the global env to get dati!!
        # get good start guesses
        hi$counts <- hi$counts[-1]
        idx       <- which(hi$counts==max(hi$counts))[1]
        more.frequent.speed <- hi$mids[idx] # assumed to be for fishing
 
        a.flag <- FALSE
        if(is.na(more.frequent.speed) | length(x)==0){ a.flag <- TRUE 
        } else {
         while(more.frequent.speed > 700 || more.frequent.speed < 100){
         hi$counts <- hi$counts[-idx]
         hi$mids <-  hi$mids[-idx]
         idx       <- which(hi$counts==max(hi$counts))[1]
         more.frequent.speed <- hi$mids[idx]
         if(is.na(more.frequent.speed)){
          #=> for very few cases we have no speed found within 100-700 e.g. for the UNK gear
          a.flag <- TRUE ; break 
          }
        }
       }  
    if(!a.flag){  
    start.lower.bound <- ifelse(more.frequent.speed-200<= min(dati$x),
                                         min(dati$x)+100, more.frequent.speed-200)
    start.upper.bound <- ifelse(more.frequent.speed+200>= max(dati$x),
                                         max(dati$x)-100, more.frequent.speed+200)
    psi    <- list(x= c(start.lower.bound, start.upper.bound) )
    psi$x[1]    <- dati$x [dati$x<=psi$x[1]] [length(dati$x [dati$x<=psi$x[1]])]   # get the bound value of the start guess from dati$x
    psi$x[2]    <- dati$x [dati$x>psi$x[2]] [1]
    o <- 1 ; class(o) <- "try-error" ; count <- 0  ; bound1 <-NULL ; bound2 <- NULL;
   while(class(o)=="try-error"){
    count <- count+1
    o <- try(
         segmented(lm(y~x, data=dati) , seg.Z=~x , psi=psi, control= seg.control(display = FALSE, it.max=50, h=1)), # with 2 starting guesses
         silent=TRUE) # the second breakpoint is quite uncertain and could lead to failure so...
   if(!"try-error" %in% class(o)) break else psi <- list(x=c(psi$x[1],psi$x[2]-20)) # searching decreasing by 20 each time
   if(count>10) {bound1 <- start.lower.bound; bound2 <- start.upper.bound ; cat(paste("failure of the segmented regression for",a.vesselid,gr,"\n...")); break}
   }
  if(is.null(bound1)) bound1 <- o$psi[order(o$psi[,"Est."])[1],"Est."] -20 # -20 hard to justify...
  if(is.null(bound2)) bound2 <- o$psi[order(o$psi[,"Est."])[2],"Est."] +20

 if(general$visual.check){
   X11()
   par(mfrow=c(2,1))
   if(!"try-error" %in% class(o)){
     plot(dati$x/100, o$fitted.values, type="l",
         ylab="cumulative distrib.", xlab="Knots",
           main=paste("segmented regression  - ",a.vesselid))
    } else{plot(0,0,type="n")}
   #plot(hi)
   #points(dati$x,dati$y)
   tmp <- as.numeric(as.character(sort(xxx$speed)))
   his <- hist(tmp, nclass=100, main=paste(general$speed,"speed between consecutive points"),
                      xlab="apparent speed [knots]", plot=TRUE)
   if(!is.null(bound1)) abline(v=bound1/100,col=2)
   if(!is.null(bound2)) abline(v=bound2/100,col=2)
   if(!is.null(bound1)) text(bound1/100, median(his$counts), signif(bound1,3), col=2, cex=2)
   if(!is.null(bound2)) text(bound2/100,  median(his$counts), signif(bound2,3), col=2, cex=2)
   # save the panel plot
   savePlot(filename = file.path(general$output.path,
      paste(unique(a.vesselid),"-detected_speed_span_for_feffort-", general$a.year,"-",gr, sep="")),
          type = c("wmf"), device = dev.cur(), restoreConsole = TRUE)
  dev.off()
  }

  # so,
  bound1 <- bound1  / 100  # transform back
  bound2 <- bound2  / 100  # transform back
  xxx$speed <- replace(xxx$speed, is.na(xxx$speed), 0) # debug 0/0
 


  # maybe you want to only keep the upper bound
  # and then set the lower one for a particular gear?
  if(gr %in% gears.to.force){
                bound1 <- force.lower.bound
              cat(paste(gr," force lower bound at", force.lower.bound, "knots\n"))
             }

  # then, assign...
  xxx[xxx$speed < bound1, "SI_STATE"]                        <- 2 # steaming
  xxx[xxx$speed >= bound1 & xxx$speed < bound2, "SI_STATE"]  <- 1 # fishing
  xxx[xxx$speed >= bound2 , "SI_STATE"]                      <- 2 # steaming
  tacsat.this.vessel[tacsat.this.vessel$LE_GEAR==gr, "SI_STATE"] <- xxx$SI_STATE # output
  tacsat.this.vessel[tacsat.this.vessel$LE_GEAR==gr,"bound1"] <- bound1
  tacsat.this.vessel[tacsat.this.vessel$LE_GEAR==gr,"bound2"] <- bound2
  cat(paste(gr," lower", general$speed,"speed bound:",round(bound1,1),"knots\n"))
  cat(paste(gr," upper", general$speed,"speed bound:",round(bound2,1),"knots\n"))
  } else{
     tacsat.this.vessel[tacsat.this.vessel$LE_GEAR==gr, "SI_STATE"] <- 2
     #=> end a.flag: in case of very few records for this gear...
     } 
  } # end gr


  # clean up by removing no longer used columns
  tacsat.this.vessel <- tacsat.this.vessel[, !colnames(tacsat.this.vessel) %in%
                         c('speed','diff.time.mins')]

  tacsat[tacsat$VE_REF==a.vesselid,] <- tacsat.this.vessel # output
  } # end of a.vesselid

 
  # plot for speed bound distribution
  tacsat2 <- tacsat[tacsat$VE_REF %in% vessels & !duplicated(tacsat[,c("VE_REF","LE_GEAR")]),]
  b1 <- tapply(tacsat2$bound1, tacsat2$LE_GEAR, mean, na.rm=TRUE)
  b2 <- tapply(tacsat2$bound2, tacsat2$LE_GEAR, mean, na.rm=TRUE)
  cat(paste("lower speed bound mean:",round(b1,1),"knots\n"))
  cat(paste("upper speed bound mean:",round(b2,1),"knots\n"))
  if(general$visual.check){
    X11()
    op <- par(no.readonly = TRUE)
    par(mfrow=c(1,2))
    boxplot(bound1 ~ LE_GEAR, data=tacsat2, ylim=c(0,8))
    boxplot(bound2 ~ LE_GEAR, data=tacsat2, ylim=c(0,8))
    par(op)
    dev.off()
  }
  # save
  write.table(tacsat2[,c("VE_REF","LE_GEAR","bound1","bound2")], append=TRUE,
   file = file.path(general$output.path,
      paste("speed_bounds_per_vessel_per_gear_", general$a.year,".txt", sep="")),
        quote=FALSE, col.names=FALSE, row.names=FALSE)


return(tacsat[tacsat$VE_REF %in% vessels, !colnames(tacsat) %in% c('bound1','bound2')])
}

