#' shortcut for as.numeric
#' 
#' Change the class of an object to numeric
#' 
#' 
#' @param x object to turn into numeric
#' @return as.numeric attempts to coerce its argument to numeric type
#' @author Niels T. Hintzen
#' @seealso \code{\link{as.numeric}}
#' @references EU Lot 2 project
#' @examples
#' 
#' as.numeric("5")   #returns the character 5 as class 'numeric'
#' an("5")           #returns the character 5 also as class 'numeric'
#' 
#' @export an
`an` <-
function(x){return(as.numeric(x))}

