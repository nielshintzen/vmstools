#' Sorting Tacsat dataframe
#' 
#' Sort the Tacsat data first by vessel, then by date, speed and heading. Needs
#' to be in this order to be effectively used in other EU lot 2 project generic
#' functions.
#' 
#' 
#' @param dat tacsat dataframe
#' @note Uses library(doBy)
#' @author Niels T. Hintzen
#' @seealso \code{\link{filterTacsat}}
#' @references EU lot 2 project
#' @examples
#' 
#' data(tacsat)
#' require(doBy)
#' 
#'   #Sort the Tacsat data
#' tacsat     <- sortTacsat(tacsat)
#' 
#' @export sortTacsat
`sortTacsat` <-
function(dat){
require(doBy)

if(!"SI_DATIM" %in% colnames(dat)) dat$SI_DATIM  <- as.POSIXct(paste(dat$SI_DATE,  dat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M")

  #Sort the tacsat data first by ship, then by date
if("VE_REF" %in% colnames(dat)) dat <- orderBy(~VE_REF+SI_DATIM,data=dat)
if("OB_REF" %in% colnames(dat)) dat <- orderBy(~OB_REF+SI_DATIM,data=dat)

return(dat)}

                                                
