formatTacsat <- function(x){
  x$VE_REF    <- ac(x$VE_REF)
  x$SI_LATI   <- an(x$SI_LATI)
  x$SI_LONG   <- an(x$SI_LONG)
  x$SI_DATE   <- ac(x$SI_DATE)
  x$SI_TIME   <- ac(x$SI_TIME)
  x$SI_SP     <- an(x$SI_SP)
  x$SI_HE     <- an(x$SI_HE)
  if("SI_HARB" %in% colnames(x))  x$SI_HARB   <- ac(x$SI_HARB)
  if("SI_STATE" %in% colnames(x)) x$SI_STATE  <- ac(x$SI_STATE)
  if("SI_FT" %in% colnames(x))    x$SI_FT     <- ac(x$SI_FT)
  #Get rid of NAs in the long and lats
  x <- x[!is.na(x$SI_LATI),]
  x <- x[!is.na(x$SI_LONG),]
  return(x)
}