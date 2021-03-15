#' Turn eflalo dataset into relational format
#' 
#' Turn the column setup of eflalo into a row setup where each species catch
#' has its own row
#' 
#' May take a long time for long eflalo datasets
#' 
#' @param x Dataframe with eflalo data and eflalo format
#' @author Niels T. Hintzen
#' @seealso \code{\link{formatEflalo}},\code{\link{readEflalo}}
#' @references EU Lot 2 project
#' @examples
#' 
#' data(eflalo)
#' eflalo    <- eflalo[1:20,]
#' eflaloRel <- eflalo2relational(eflalo)
#' 
#' @export eflalo2relational
eflalo2relational <- function(x){
  idxkg             <- grep("KG",colnames(x))
  idxeuro           <- grep("EURO",colnames(x))
  selectColsKG      <- apply(x[,idxkg],  1,FUN=function(y){idx <- which(y >0);
                                                           return(y[idx])})
  selectColsEURO    <- apply(x[,idxeuro],1,FUN=function(y){idx <- which(y >0);
                                                           return(y[idx])})


  idxNonkgeur       <- which(!1:length(colnames(x)) %in% kgeur(colnames(x)))
  xRelationalList   <- lapply(as.list(1:length(selectColsKG)),function(y){
                               newX              <- data.frame(x[y,idxNonkgeur])
                               newX$LE_SP        <- strsplit(names(selectColsKG[[y]][1]),"_")[[1]][3]
                               newX$LE_KG        <- selectColsKG[[y]][1]
                               newX$LE_EURO      <- selectColsEURO[[y]][1]
                               if(length(selectColsKG[[y]])>1){
                                 for(j in 2:length(selectColsKG[[y]])){
                                   newX <- rbind(newX,
                                                 data.frame(x[y,idxNonkgeur],LE_SP=strsplit(names(selectColsKG[[y]][j]),"_")[[1]][3],
                                                                             LE_KG=selectColsKG[[y]][j],
                                                                             LE_EURO=selectColsEURO[[y]][j]))
                                 }
                               }
                               return(newX)})
  xRelational       <- do.call(rbind,xRelationalList)
return(xRelational)}
#eflalo2 <- eflalo2relational(eflalo)

