#' Produce quoted list for use within data.table
#' 
#' Produce the quoted list format from a character vector which is used as
#' input for the data.table format
#' 
#' character vector could contain e.g. colnames(x) or c("COL1","COL2","COL3")
#' where COL1 etc are column names
#' 
#' @param \dots character vector
#' @return Returnes quoted list
#' @note With great thanks to original author.
#' @author Niels T. Hintzen
#' @seealso See \code{\link{data.table}} documentation
#' @references See \code{\link{data.table}} documentation
#' @examples
#' 
#' c.listquote(c("COL1","COL2","COL3")) #Returns: list(COL1,COL2,COL3)
#' 
#' @export c.listquote
c.listquote <- function( ... ) {

   args <- as.list( match.call()[ -1 ] )
   lstquote <- list( as.symbol( "list" ) );
   for ( i in args ) {
      # Evaluate expression in parent eviron to see what it refers to
      if ( class( i ) == "name" || ( class( i ) == "call" && i[[1]] != "list" ) ) {
         i <- eval( substitute( i ), sys.frame( sys.parent() ) )
      }
      if ( class( i ) == "call" && i[[1]] == "list" ) {
         lstquote <- c( lstquote, as.list( i )[ -1 ] )
      }
      else if ( class( i ) == "character" )
      {
         for ( chr in i ) {
            lstquote <- c( lstquote, list( parse( text=chr )[[1]] ) )
         }
      }
      else
         stop( paste( "[", deparse( substitute( i ) ), "] Unknown class [", class( i ), "] or is not a list()", sep="" ) )
   }
   return( as.call( lstquote ) )
}
