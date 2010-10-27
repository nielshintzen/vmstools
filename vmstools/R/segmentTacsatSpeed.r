# do the distinction between fishing and non-fishing
# by automatic detection of the fishing peak
# required: library 'segmented'

segmentSpeedTacsat <- function(tacsat, vessels="DNK100"){

  library(segmented)
  tacsat[,"bound1"] <- NA
  tacsat[,"bound2"] <- NA

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



  for(a.vesselid in vessels){

  tacsat.this.vessel <- tacsat[tacsat$VE_REF %in% a.vesselid, ]
  
  Sys.setlocale("LC_TIME", "english")
  ctime <- strptime(  paste(tacsat.this.vessel$SI_DATE, tacsat.this.vessel$SI_TIME) ,
                                   "%e/%m/%Y %H:%M" )
  tacsat.this.vessel <- cbind.data.frame(tacsat.this.vessel, date.in.R=ctime)

  diff.time <- tacsat.this.vessel[-nrow(tacsat.this.vessel),"date.in.R"] -
                       tacsat.this.vessel[-1,"date.in.R"]
  tacsat.this.vessel$diff.time.mins <- c(0, as.numeric(diff.time, units="mins"))

  # add a apparent speed colunm (nautical miles per hour)
  tacsat.this.vessel$apparent.speed <- signif(
            c(0, distAB.f(A= tacsat.this.vessel[-nrow(tacsat.this.vessel),c("SI_LATI","SI_LONG")],
              B= tacsat.this.vessel[-1,c("SI_LATI","SI_LONG")], .unit="nm") / (abs(xx[-1,]$diff.time.mins)/60)),
              3)

  # cleaning irrealistic points
  tacsat.this.vessel$apparent.speed <-
     replace(tacsat.this.vessel$apparent.speed, is.na(tacsat.this.vessel$apparent.speed), 0)
  tacsat.this.vessel <- tacsat.this.vessel[tacsat.this.vessel$apparent.speed < 30,]

  for (gr in levels(tacsat$LE_GEAR)){
    xxx <- tacsat.this.vessel$LE_GEAR[tacsat.this.vessel$LE_GEAR==gr,] # input

    x <- as.numeric(as.character(sort(xxx$apparent.speed))) *100   # multiply by factor 100 because integer needed
    hi <- hist(x, nclass=30,plot=FALSE) # remove irealistic speed at the mean time

    y <- c(1:length(sort(xxx$apparent.speed))) # sort in increasing order
    y <- y[x>100 & x<1000] # just removing the 0, and the larger speeds we 100% know it is steaming
    x <- x[x>100 & x<1000]
    dati   <- data.frame(x=x,y=y)
    dati$x <- as.integer(dati$x) # integer needed
    psi    <- list(x= quantile(dati$x,probs=c(0.05,0.5))  )
    assign('dati', dati, env=.GlobalEnv) # DEBUG segmented()...pfff! this function looks in the global env to get dati!!
    # get good start guesses
    hi$counts <- hi$counts[-1]
    idx       <- which(hi$counts==max(hi$counts))[1]
    more.frequent.speed <- hi$mids[idx] # assumed to be for fishnig
    while(more.frequent.speed > 700 || more.frequent.speed < 100){
      hi$counts <- hi$counts[-idx]
      hi$mids <-  hi$mids[-idx]
      idx       <- which(hi$counts==max(hi$counts))[1]
      more.frequent.speed <- hi$mids[idx]
      }
    start.lower.bound <- ifelse(more.frequent.speed-200<= min(dati$x), min(dati$x)+100, more.frequent.speed-200)
    start.upper.bound <- ifelse(more.frequent.speed+200>= max(dati$x), max(dati$x)-100, more.frequent.speed+200)
    psi    <- list(x= c(start.lower.bound, start.upper.bound) )
    psi$x[1]    <- dati$x [dati$x<=psi$x[1]] [length(dati$x [dati$x<=psi$x[1]])]   # get the bound value of the start guess from dati$x
    psi$x[2]    <- dati$x [dati$x>psi$x[2]] [1]
    o <- 1 ; class(o) <- "try-error" ; count <- 0  ; bound1 <-NULL ; bound2 <- NULL;
   while(class(o)=="try-error"){
    count <- count+1
    o <- try(
         segmented(lm(y~x, data=dati) , seg.Z=~x , psi=psi, control= seg.control(display = FALSE, it.max=50, h=1)), # with 2 starting guesses
         silent=TRUE) # the second breakpoint is quite uncertain and could lead to failure so...
   if(class(o)!="try-error") break else psi <- list(x=c(psi$x[1],psi$x[2]-20)) # searching decreasing by 100 each time
   if(count>10) {bound1 <- start.lower.bound; bound2 <- start.upper.bound ; cat(paste("failure of the segmented regression for",a.vesselid,gr,"\n...")); break}
   }
  if(is.null(bound1)) bound1 <- o$psi[order(o$psi[,"Est."])[1],"Est."] -20 # -20 hard to justify...
  if(is.null(bound2)) bound2 <- o$psi[order(o$psi[,"Est."])[2],"Est."] +20

 if(TRUE){
   windows()
   par(mfrow=c(3,1))
   if(class(o)!="try-error"){
     plot(dati$x/100,o$fitted.values, type="l",ylab="cumulative distrib.", xlab="Knots", main=paste("segmented regression  - ",a.vesselid))
    }
   plot(hi)
   #points(dati$x,dati$y)
   tmp <- as.numeric(as.character(sort(xxx$apparent.speed)))
   hist(tmp, nclass=100, main="apparent speed between consecutive points", xlab="apparent speed [knots]")
   if(!is.null(bound1)) abline(v=bound1/100,col=2)
   if(!is.null(bound2)) abline(v=bound2/100,col=2)
   if(!is.null(bound1)) text(bound1/100, 1, signif(bound1,3), col=2)
   if(!is.null(bound2)) text(bound2/100, 1, signif(bound2,3), col=2)
   # save the panel plot
   savePlot(filename = file.path(general$output.path,
      paste(unique(a.vesselid),"-detected_speed_span_for_feffort-", general$a.year,"-",a.vesselid,"-",gr, sep="")),
          type = c("wmf"), device = dev.cur(), restoreConsole = TRUE)
  dev.off()
  }

  # so,
  bound1 <- bound1  / 100  # re-transform
  bound2 <- bound2  / 100   # re-transform
  xxx$apparent.speed <- replace(xxx$apparent.speed, is.na(xxx$apparent.speed), 0) # debug 0/0

  if(gr %in% c('GNS','GND','GNC','GNF',
             'GTR','GTN','GEN','GN')) bound1 <- 0.5 # special case for gillnetters= no lower bound
  xxx[xxx$apparent.speed < bound1, "state"]                       <- 2 # steaming
  xxx[xxx$apparent.speed >= bound1 & xxx$apparent.speed < bound2, "state"]  <- 1 # fishing
  xxx[xxx$apparent.speed >= bound2 , "state"]                     <- 2 # steaming
  xxx[,"bound1"] <- bound1
  xxx[,"bound2"] <- bound2
  tacsat.this.vessel[tacsat.this.vessel$gear==gr,] <- xxx # output
  } # end gr

  } # end of a.vesselid

return(tacsat.this.vessel)
}

