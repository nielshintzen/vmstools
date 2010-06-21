ac <- function(x){return(as.character(x))}
an <- function(x){return(as.numeric(x))}
af <- function(x){return(as.factor(x))}

  #Get the number of digits
getndp <- function(x, tol=2*.Machine$double.eps)
{
  ndp <- 0
  while(!isTRUE(all.equal(x, round(x, ndp), tol=tol))) ndp <- ndp+1
  if(ndp > -log10(tol)) warning("Tolerance reached, ndp possibly
underestimated.")
  ndp
}

  #Calculate the 2-dimensional normal distribution
N1p0 <- function(x,mu,sig,p){
          aa<- c((1/(sig*sqrt(2*pi)))*exp(-((x-(mu+p))^2)/(2*sig^2)))
          return(aa)
      }
      
N1 <- function(x,sig){
          aa<- c((1/(sig*sqrt(2*pi)))*exp(-((x)^2)/(2*sig^2)))
          return(aa)
      }