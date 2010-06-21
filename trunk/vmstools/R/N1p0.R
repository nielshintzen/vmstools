`N1p0` <-
function(x,mu,sig,p){
          aa<- c((1/(sig*sqrt(2*pi)))*exp(-((x-(mu+p))^2)/(2*sig^2)))
          return(aa)
      }

