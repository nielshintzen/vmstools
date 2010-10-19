formatEflalo2 <- function(x){
  x$VE_REF        <- ac(x$VE_REF)
  x$VE_FLT        <- ac(x$VE_FLT)
  x$VE_COU        <- ac(x$VE_COU)
  x$VE_LEN        <- an(ac(x$VE_LEN))
  x$VE_KW         <- an(ac(x$VE_KW))
  x$VE_TON        <- an(ac(x$VE_TON))
  x$FT_REF        <- ac(x$FT_REF)
  x$FT_DCOU       <- ac(x$FT_DCOU)
  x$FT_DHAR       <- ac(x$FT_DHAR)
  x$FT_DDAT       <- ac(x$FT_DDAT)
  x$FT_DTIME      <- ac(x$FT_DTIME)
  x$FT_LCOU       <- ac(x$FT_LCOU)
  x$FT_LHAR       <- ac(x$FT_LHAR)
  x$FT_LDAT       <- ac(x$FT_LDAT)
  x$FT_LTIME      <- ac(x$FT_LTIME)
  x$LE_CDAT       <- ac(x$LE_CDAT)
  x$LE_STIME      <- ac(x$LE_STIME)
  x$LE_ETIME      <- ac(x$LE_ETIME)
  x$LE_SEQNUM     <- an(ac(x$LE_SEQNUM))
  x$LE_GEAR       <- ac(x$LE_GEAR)
  x$LE_MSZ        <- an(ac(x$LE_MSZ))
  x$LE_RECT       <- ac(x$LE_RECT)
  x$LE_MET_level6 <- ac(x$LE_MET_level6)
  x$LE_UNIT       <- ac(x$LE_UNIT)
  x$LE_EFF        <- an(ac(x$LE_EFF))
  x$LE_EFF_VMS    <- an(ac(x$LE_EFF_VMS))
  for(i in c(grep("_KG_",colnames(x)),grep("_EURO_",colnames(x)))) x[,i] <- an(ac(x[,i]))
  return(x)
}




