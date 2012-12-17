
interCubicHermiteSpline <- function(spltx,spltCon,res,params,headingAdjustment){

      #Formula of Cubic Hermite Spline
      t   <- seq(0,1,length.out=res)
      F00 <- 2*t^3 -3*t^2 + 1
      F10 <- t^3-2*t^2+t
      F01 <- -2*t^3+3*t^2
      F11 <- t^3-t^2

      #Making tacsat dataset ready
      spltx[spltCon[,1],"SI_HE"][which(is.na(spltx[spltCon[,1],"SI_HE"]))] <- 0
      spltx[spltCon[,2],"SI_HE"][which(is.na(spltx[spltCon[,2],"SI_HE"]))] <- 0

      #Heading at begin point in degrees
      Hx0 <- sin(spltx[spltCon[,1],"SI_HE"]/(180/pi))
      Hy0 <- cos(spltx[spltCon[,1],"SI_HE"]/(180/pi))

      #Heading at end point in degrees
      Hx1 <- sin(spltx[spltCon[,2]-headingAdjustment,"SI_HE"]/(180/pi))
      Hy1 <- cos(spltx[spltCon[,2]-headingAdjustment,"SI_HE"]/(180/pi))

      #Start and end positions
      Mx0 <- spltx[spltCon[,1],"SI_LONG"]
      Mx1 <- spltx[spltCon[,2],"SI_LONG"]
      My0 <- spltx[spltCon[,1],"SI_LATI"]
      My1 <- spltx[spltCon[,2],"SI_LATI"]

      #Corrected for longitude lattitude effect
      Hx0 <- Hx0 * params$fm * spltx[spltCon[,1],"SI_SP"] /((params$st[2]-params$st[1])/2+params$st[1])
      Hx1 <- Hx1 * params$fm * spltx[spltCon[,2],"SI_SP"] /((params$st[2]-params$st[1])/2+params$st[1])
      Hy0 <- Hy0 * params$fm * lonLatRatio(spltx[spltCon[,1],"SI_LONG"],spltx[spltCon[,1],"SI_LATI"]) * spltx[spltCon[,1],"SI_SP"]/((params$st[2]-params$st[1])/2+params$st[1])
      Hy1 <- Hy1 * params$fm * lonLatRatio(spltx[spltCon[,2],"SI_LONG"],spltx[spltCon[,2],"SI_LATI"]) * spltx[spltCon[,2],"SI_SP"]/((params$st[2]-params$st[1])/2+params$st[1])

      #Get the interpolation
      fx  <- outer(F00,Mx0,"*")+outer(F10,Hx0,"*")+outer(F01,Mx1,"*")+outer(F11,Hx1,"*")
      fy  <- outer(F00,My0,"*")+outer(F10,Hy0,"*")+outer(F01,My1,"*")+outer(F11,Hy1,"*")

      #Create output format
      intsx   <- lapply(as.list(1:nrow(spltCon)),function(x){
                    matrix(rbind(spltx$ID[spltCon[x,]],cbind(fx[,x],fy[,x])),ncol=2,
                           dimnames=list(c("startendVMS",seq(1,res,1)),c("x","y")))})
  return(intsx)}

