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

