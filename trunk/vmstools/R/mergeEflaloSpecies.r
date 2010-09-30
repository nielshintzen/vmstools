 # reduce the size of the eflalo data by merging species
 # threshold in euros
mergeEflaloSpecies <- function(eflalo, threshold.value=3000000){
    get.sp <- function (nm) unlist(lapply(strsplit(nm, split="_"), function(x) x[3]))
    nm         <- names(eflalo)
    idx.col.w  <- grep('KG', nm) # index columns weight
    idx.col.v  <- grep('EURO', nm)  # index columns value
    v          <- apply(eflalo [,idx.col.v], 2, sum, na.rm=TRUE)
    sp.to.merge.v <- names(v)[v < threshold.value] # e.g. threshold in euros
    sp.to.merge.w <- paste('LE_KG_', get.sp(sp.to.merge.v ), sep='')
    sp.to.merge.w <- sp.to.merge.w[sp.to.merge.w != 'LE_KG_MZZ']
    sp.to.merge.v <- sp.to.merge.v[sp.to.merge.v != 'LE_EURO_MZZ']
    OTH.w    <- apply(eflalo[,sp.to.merge.w], 1, sum, na.rm=TRUE)
    OTH.v    <- apply(eflalo[,sp.to.merge.v], 1, sum, na.rm=TRUE)
    eflalo  <- eflalo [, !nm %in% c(sp.to.merge.w, sp.to.merge.v)]
    eflalo[,'LE_KG_MZZ'] <- replace(eflalo[,'LE_KG_MZZ'], is.na(eflalo[,'LE_KG_MZZ']), 0)
    eflalo[,'LE_KG_MZZ']   <- eflalo[,'LE_KG_MZZ'] + OTH.w
    eflalo[,'LE_EURO_MZZ'] <- replace(eflalo[,'LE_EURO_MZZ'], is.na(eflalo[,'LE_EURO_MZZ']), 0)
    eflalo[,'LE_EURO_MZZ'] <- eflalo[,'LE_EURO_MZZ'] + OTH.v
  return(eflalo)
  }
              