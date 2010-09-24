 # reduce the size of the eflalo data by merging species
 # threshold in kilo
mergeEflaloSpecies <- function(eflalo, threshold=100000){
    get.sp <- function (nm) unlist(lapply(strsplit(nm, split="_"), function(x) x[3]))
    nm         <- names(eflalo)
    idx.col.w  <- grep('KG', nm) # index columns weight
    idx.col.v  <- grep('EURO', nm)  # index columns value
    w          <- apply(eflalo [,idx.col.w], 2, sum, na.rm=TRUE)
    sp.to.merge.w <- names(w)[w < threshold] # e.g. threshold: 100 tons inside this year
    sp.to.merge.v <- paste('LE_EURO_', get.sp(sp.to.merge.w ), sep='')
    sp.to.merge.w <- sp.to.merge.w[sp.to.merge.w != 'LE_KG_MZZ']
    sp.to.merge.v <- sp.to.merge.v[sp.to.merge.v != 'LE_EURO_MZZ']
    OTH.w    <- apply(eflalo2[,sp.to.merge.w], 1, sum, na.rm=TRUE)
    OTH.v    <- apply(eflalo2[,sp.to.merge.v], 1, sum, na.rm=TRUE)
    eflalo  <- eflalo [, !nm %in% c(sp.to.merge.w, sp.to.merge.v)]
    eflalo[,'LE_KG_MZZ']   <- eflalo2[,'LE_KG_MZZ'] + OTH.w
    eflalo[,'LE_EURO_MZZ'] <- eflalo2[,'LE_EURO_MZZ'] + OTH.v
  return(eflalo)
  }
