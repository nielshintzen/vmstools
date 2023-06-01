#' Get Number of Decimal Places
#' 
#' Return the number of decimal places of a 'numeric'.
#' 
#' 
#' @param x Number to find decimal places from
#' @param tol Tolerance to use
#' @note Function not created under EU lot 2 project but found under R-help.
#' Please look there for credits.
#' @author See R-help pages
#' @references EU lot 2 project
#' @examples
#' 
#' getndp(5.677) #result: 3
#' 
#' @export getndp
`getndp` <-
function(x, tol=2*.Machine$double.eps)
{
  ndp <- 0
  while(!isTRUE(all.equal(x, round(x, ndp), tol=tol))) ndp <- ndp+1
  if(ndp > -log10(tol)) warning("Tolerance reached, ndp possibly
underestimated.")
  ndp
}

