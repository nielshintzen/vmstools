estimateAggregation <- function(x,trustZeros=T){

    if(all(names(x) %in% c(0,1)))
      stop("No aggregation in binomial data")
    if(trustZeros==T & !0 %in% an(names(x)))
      stop("trustZero cannot be TRUE when no zeros are given")
    if(trustZeros==F & 0 %in% an(names(x)))
      x       <- x[-which(names(x)==0)]
    startMu   <- sum(an(names(x))*x)/sum(x)
    startSize <- startMu/0.2
    dummy     <- 1e-8
    if(0 %in% an(names(x)))
      funNB   <- function(par,x){
                 return(sum(-1*
                            dnbinom(an(names(x)),size=(abs(par[1])+dummy)/abs(par[2]),mu=abs(par[1])+dummy,log=T) * x,
                        na.rm=T))}
    if(!0 %in% an(names(x)))
      funNB   <- function(par,x){
                 return(sum(-1*
                            (dnbinom(an(names(x)),size=(abs(par[1])+dummy)/abs(par[2]),mu=abs(par[1])+dummy,log=T) -
                            log(1-dnbinom(0,size=(abs(par[1])+dummy)/abs(par[2]),mu=abs(par[1])+dummy,log=F))) * x,
                        na.rm=T))}
    ret       <- abs(optim(par=c(startMu,startSize),x=x,fn=funNB)$par)
    names(ret)<- c("Mean trawling intensity","Aggregation")
    return(ret)}
    
estimateZeros <- function(x){
                    mubeta      <- estimateAggregation(x)
                    zeros       <- sum(an(names(x))*x)/mubeta[1] - sum(x)
                 return(zeros)}
                 
predictAggDistri <- function(x,trustZeros=T){
                        mubeta      <- estimateAggregation(x,trustZeros)
                        if(trustZeros == F){
                          zeros     <- sum(an(names(x))*x)/mubeta[1] - sum(x)
                        } else {
                          zeros     <- numeric()
                        }
                        if(0 %in% an(names(x)))
                          vals      <- an(names(x))
                        if(!0 %in% an(names(x)))
                          vals      <- c(0,an(names(x)))

                        pred        <- dnbinom(vals,size=mubeta[1]/mubeta[2],mu=mubeta[1])*sum(c(zeros,x))
                        ret         <- cbind(Observed = round(c(zeros,x)),Predicted = round(pred))
                        dimnames(ret)[[1]] <- vals
                    return(ret)}

