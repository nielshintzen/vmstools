#' Row bind two tacsat sets together
#' 
#' Row binds two tacsat sets together while taking differences in column names
#' into account
#' 
#' 
#' @param set1 Tacsat set 1
#' @param set2 Tacsat set 2
#' @author Niels T. Hintzen
#' @seealso \code{\link{rbindEflalo}}, \code{\link{do.call}}
#' @references EU Lot 2 project
#' @examples
#' 
#' data(tacsat)
#' set1 <- tacsat
#' set2 <- tacsat[seq(1,100,5),]
#' 
#' combined <- rbindTacsat(set1,set2)
#' 
#' 
#' @export rbindTacsat
rbindTacsat <- function(set1,set2){
  cln1  <- colnames(set1)
  cln2  <- colnames(set2)
  if(any(duplicated(cln1)==TRUE) || any(duplicated(cln2)==TRUE)) stop("Duplicate column names in datasets")
  idx1  <- which(is.na(pmatch(cln1,cln2))==TRUE)
  idx2  <- which(is.na(pmatch(cln2,cln1))==TRUE)
  
  if(length(idx1)>0){
    for(i in idx1) set2 <- cbind(set2,NA)
    colnames(set2) <- c(cln2,cln1[idx1])}
  if(length(idx2)>0){
    for(i in idx2) set1 <- cbind(set1,NA)
    colnames(set1) <- c(cln1,cln2[idx2])}
  cln1  <- colnames(set1)
  cln2  <- colnames(set2)
  mtch  <- pmatch(cln1,cln2)
  if(any(is.na(mtch))==TRUE) stop("Cannot find nor create all matching column names")
  set3  <- rbind(set1,set2[,cln2[mtch]])
return(set3)}
  



