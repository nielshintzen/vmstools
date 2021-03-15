#' shortcut for as.numeric(as.character())
#' 
#' Change the class of an object from factor to numeric
#' 
#' 
#' @param x object to turn from factor into numeric
#' @return as.numeric attempts to coerce its argument to numeric type
#' @author Francois Bastardie
#' @seealso \code{\link{as.numeric}}, \code{\link{as.character}}
#' @references EU Lot 2 project
#' @examples
#' 
#' 
#' res <- as.factor(5.1)
#' an(res)   #returns 1
#' anf(res)  #returns the original 5.1
#' 
#' @export anf
  'anf' <-
  function(x) as.numeric(as.character(x)) # alias to convert factors
