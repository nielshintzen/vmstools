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
