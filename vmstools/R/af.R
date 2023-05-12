#' shortcut for as.factor
#' 
#' Change the class of an object to factor
#' 
#' 
#' @param x object to turn into factor
#' @return as.factor attempts to coerce its argument to factor type
#' @author Niels T. Hintzen
#' @seealso \code{\link{as.factor}}
#' @references EU Lot 2 project
#' @examples
#' 
#' as.factor(5)    #returns the number 5 as class 'factor'
#' af(5)           #returns the number 5 also as class 'factor'
#' 
#' @export af
`af` <-
function(x){return(as.factor(x))}

#hello world
