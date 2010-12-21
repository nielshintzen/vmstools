   # add a year.quarter column
   addQuarter <- function(obj){
      if(!"SI_DATE" %in% colnames(obj))
          stop("a SI_DATE column is required")
     obj$month <- factor(format(as.POSIXct(obj$SI_DATE), "%m"))  # init
     obj$quarter <- obj$month # init
     levels(obj$quarter) <- c(1,1,1,2,2,2,3,3,3,4,4,4)
   return(obj)
   }
   # all.merged <- addQuarter(all.merged)
