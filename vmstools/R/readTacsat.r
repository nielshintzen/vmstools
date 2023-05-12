#' Read Tacsat data into R
#' 
#' Reading tacsat data from delimited file into R, checking for obligatory
#' columns and formatting all data columns.
#' 
#' 
#' @param file file path + file name
#' @param sep delimiter used in file (default to ',')
#' @param dec decimal notation used in file (default to '.')
#' @return Returns the formatted Tacsat dataset
#' @author Niels T. Hintzen
#' @seealso \code{readEflalo()}
#' @references EU lot 2 project
#' @examples
#' 
#' data(tacsat)
#' dir.create("C:/tmpTacsat")
#' 
#' #temporarily write tacsat file to disk to thereafter read it back in again
#' tacsat$SI_HE[which(tacsat$SI_HE>360 | tacsat$SI_HE < 0)] <- NA
#' tacsat$SI_LONG[which(tacsat$SI_LONG < -180 | tacsat$SI_LONG > 180)] <- NA
#' tacsat$SI_LATI[which(tacsat$SI_LATI < -90 | tacsat$SI_LATI > 90)] <- NA
#' write.table(tacsat, file = "C:/tmpTacsat/tacsat.csv", quote = TRUE, sep = ",",
#'         eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE)
#' 
#' #Read in tacsat file
#' tacsat <- readTacsat("C:/tmpTacsat/tacsat.csv")
#' 
#' 
#' 
#' @export readTacsat
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
