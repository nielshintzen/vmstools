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