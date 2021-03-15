#' shortcut for as.character
#' 
#' Change the class of an object to character
#' 
#' 
#' @param x object to turn into character
#' @return as.character attempts to coerce its argument to character type
#' @author Niels T. Hintzen
#' @seealso \code{\link{as.character}}
#' @references EU Lot 2 project
#' @examples
#' 
#' as.character(5) #returns the number 5 as class 'character'
#' ac(5)           #returns the number 5 also as class 'character'
#' 
#' @export ac
`ac` <-
function(x){return(as.character(x))}

               # Hello world !
