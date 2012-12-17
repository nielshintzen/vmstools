
interStraightLine <- function(spltx,spltCon,res){

     fx <- mapply(seq,spltx[spltCon[,1],"SI_LONG"],spltx[spltCon[,2],"SI_LONG"],length.out=res)
     fy <- mapply(seq,spltx[spltCon[,1],"SI_LATI"],spltx[spltCon[,2],"SI_LATI"],length.out=res)

     #Create output format
     intsx   <- lapply(as.list(1:nrow(spltCon)),function(x){
                   matrix(rbind(spltx$ID[spltCon[x,]],cbind(fx[,x],fy[,x])),ncol=2,
                          dimnames=list(c("startendVMS",seq(1,res,1)),c("x","y")))})
  return(intsx)}

