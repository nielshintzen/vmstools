#- Script to tune the cHs method

```


#- Load a high resolution VMS dataset
data(VMShf)
colnames(VMShf)   <- c("VE_REF","SI_LATI","SI_LONG","SI_SP","SI_HE","SI_DATIM")
VMShf$SI_DATE     <- ac(as.Date(VMShf$SI_DATIM))
VMShf$SI_TIME     <- ac(format(VMShf$SI_DATIM,format="%H:%M"))

#- Calculate speed
VMShf             <- calculateSpeed(VMShf,level="vessel",fill.na=T)

#- Generate low resolution VMS dataset from high resolution set
VMShf             <- intervalTacsat(VMShf,level="vessel",fill.na=T)
VMShf$INTVCUM     <- cumsum(VMShf$INTV)
VMSlf             <- VMShf[which(VMShf$INTVCUM %% 120 == 0),]
VMShf             <- VMShf[which(VMShf$INTVCUM %% 120 == 0)[1]:nrow(VMShf),]

#- Filter out only fishing pings
VMShf             <- subset(VMShf,SI_SPCA >= 4 & SI_SPCA <= 8.5)
VMSlf             <- subset(VMSlf,SI_SPCA >= 4 & SI_SPCA <= 8.5)

#- Define the tuning method
tunecHsInterpolation <- function(lowResTacsat,
hiResTacsat,
interval,
margin,
res,
headingAdjustment=0,
st=c(4,8.5)){

#- Define the optimization routine inside the function
optimInterpolation <- function(x){
ints <- interpolateTacsat(hiResTacsat,
interval=interval,
margin=margin,
res=res,
method="cHs",
params=list(fm=x,distscale=NA,signline=NA,st=st),
headingAdjustment=headingAdjustment)

diffs           <- diffInter(ints,hiResTacsat)[,"mean"]
pointsToRemove  <- which(diffs <= quantile(diffs,probs=c(0.01),na.rm=T) | diffs >= quantile(diffs,probs=c(0.99),na.rm=T))
meddiffs        <- median(diffs[-pointsToRemove],na.rm=T)

return(meddiffs)}

res <- optimize(optimInterpolation,interval=c(0,1),tol=0.0005)$minimum
return(res)}


#- Perform the tuning
#- There is no sensitivity to the 'st' parameter settings, as long as you use the same set of
#   two numbers in all occasions
tunecHsInterpolation(lowResTacsat=VMSlf,hiResTacsat=VMShf,interval=120,margin=3,res=100,headingAdjustment=1,st=c(4,8.5))

#- Plot the result (remember, this is an example dataset which is extremely short :-))
ints <- interpolateTacsat(VMShf,
interval=120,
margin=3,
res=100,
method="cHs",
params=list(fm=0.236955,distscale=NA,signline=NA,st=c(4,8.5)),
headingAdjustment=0)
intsTC<- interpolation2Tacsat(ints,VMSlf)
plot(intsTC$SI_LONG,intsTC$SI_LATI,type="l",xlab="Longitude",ylab="Latitude")
lines(VMShf$SI_LONG,VMShf$SI_LATI,col=2)
points(VMSlf$SI_LONG,VMSlf$SI_LATI,pch=19,cex=1,col="blue")
```