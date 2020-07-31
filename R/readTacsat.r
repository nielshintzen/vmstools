readTacsat <- function(file,sep=",",dec="."){

                    #Read the data
                  res <- read.table(file, header = TRUE,sep,dec = ".",stringsAsFactors = FALSE)

                    #Perform checks
                  if(any(!c("VE_REF","SI_LATI","SI_LONG","SI_DATE","SI_TIME","SI_SP","SI_HE") %in% colnames(res))) stop(paste(file,"needs correct header including",paste("VE_REF","SI_LATI","SI_LONG","SI_DATE","SI_TIME","SI_SP","SI_HE")))
                  if(any(res$SI_LATI > 90 | res$SI_LATI < -90,na.rm=TRUE) | any(res$SI_LONG > 180 | res$SI_LONG < -180,na.rm=TRUE)) stop("Longitudes or latitudes are out of range")
                  if(any(res$SI_HE > 360 | res$SI_HE < 0,na.rm=TRUE)) stop("Heading out of range, must be between 0 - 360")

                    #Reformat the data
                  res <- formatTacsat(res)

              return(res)}