#' Return kg and euro column index of eflalo dataset
#' 
#' Returns the index of the columns with kg and euro information from a given
#' eflalo dataset
#' 
#' 
#' @param x Colnames of eflalo dataset (or any other dataset with column names
#' as LE_KG_ and LE_EURO_)
#' @author Niels T. Hintzen
#' @references EU Lot 2 project
#' @examples
#' 
#' data(eflalo)
#' kgeur(colnames(eflalo))
#' 
#' @export kgeur
kgeur <- function(x){return(c(grep("KG",x),grep("EURO",x)))}
