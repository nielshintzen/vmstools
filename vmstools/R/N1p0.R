#' 2-dimensional Gaussian distribution with two-peak possibility
#' 
#' Computes the 2-dimensional Gaussian distribution depending on the distance
#' from the mean, and standard deviation with adjustable mean and shift from
#' mean into two peaks instead of 1
#' 
#' 
#' @param x Observed value, deviance from mean
#' @param mu Mean value
#' @param sig Standard deviation of value
#' @param p Shift from mean into two peaks
#' @note Function is called inside plotCIinterpolation()
#' @author Niels T. Hintzen
#' @references EU lot 2 project
#' @examples
#' 
#' 
#' x   <- matrix(c(5,4,3,4,5,4,3,2,3,4,3,2,1,2,3,
#'                 3,2,1,2,3,4,3,2,3,4,5,4,3,4,5),nrow=6,ncol=5,byrow=TRUE)
#' mu  <- 0
#' sig <- x^0.5
#' p   <- 0
#' 
#' res <- matrix(N1p0(x,mu,sig,p),nrow=6,ncol=5)
#' 
#'   #Plot the results
#' image(res)
#' 
#' @export N1p0
`N1p0` <-
function(x,mu,sig,p){
          aa<- c((1/(sig*sqrt(2*pi)))*exp(-((x-(mu+p))^2)/(2*sig^2)))
          return(aa)
      }

