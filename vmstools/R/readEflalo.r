#' Read Eflalo data into R
#' 
#' Reading eflalo data from delimited file into R, checking for obligatory
#' columns and formatting all data columns.
#' 
#' 
#' @param file file path + file name
#' @param sep delimiter used in file (default to ',')
#' @param dec decimal notation used in file (default to '.')
#' @return Returns the formatted Eflalo dataset
#' @author Niels T. Hintzen
#' @seealso \code{\link{readTacsat}}
#' @references EU lot 2 project
#' @examples
#' 
#' data(eflalo)
#' dir.create("C:/tmpEflalo")
#' 
#' #temporarily write tacsat file to disk to thereafter read it back in again
#' write.table(eflalo, file = "C:/tmpEflalo/eflalo.csv", quote = TRUE, sep = ",",
#'     eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = TRUE)
#' 
#' #Read in tacsat file
#' eflalo <- readEflalo("C:/tmpEflalo/eflalo.csv")
#' 
#' 
#' @export readEflalo
readEflalo <- function(file,sep=",",dec="."){

                    #Read the data
                  res <- read.table(file, header = TRUE,sep,dec = ".",stringsAsFactors = FALSE)

                    #Perform checks
                  if(any(!c("VE_REF","VE_FLT","VE_COU","VE_LEN","VE_KW","VE_TON","FT_REF","FT_DCOU","FT_DHAR","FT_DDAT","FT_DTIME",
                            "FT_LCOU","FT_LHAR","FT_LDAT","FT_LTIME","LE_ID","LE_CDAT","LE_STIME","LE_ETIME","LE_SLAT","LE_SLON",
                            "LE_ELON","LE_ELON","LE_GEAR","LE_MSZ") %in% colnames(res)))
                              stop(paste(file,"needs correct header including",paste(c("VE_REF","VE_FLT","VE_COU","VE_LEN","VE_KW","VE_TON","FT_REF","FT_DCOU","FT_DHAR","FT_DDAT","FT_DTIME",
                                  "FT_LCOU","FT_LHAR","FT_LDAT","FT_LTIME","LE_ID","LE_CDAT","LE_STIME","LE_ETIME","LE_SLAT","LE_SLON",
                                  "LE_ELON","LE_ELON","LE_GEAR","LE_MSZ"))))
                  if(length(grep("KG",colnames(res))) < 1) stop("Units of landings are not in KG")
                  if(length(grep("LE_KG",colnames(res))) < 1) stop("No landing data provided")
                  if(length(grep("LE_EURO",colnames(res))) < 1) warning("Currency used is different than 'EURO'")

                    #Reformat the data
                  res <- formatEflalo(res)

              return(res)}
