 # reduce the size of the eflalo data by merging species
 # threshold in euros
poolEflaloSpecies <- function(eflalo, threshold.value=3000000, code="MZZ"){
    if(!paste("LE_KG_", code, sep='') %in% colnames(eflalo)) {
              eflalo [, paste("LE_KG_", code, sep='')] <- NA
               eflalo [, paste("LE_EURO_", code, sep='')] <- NA
              }
    get.sp <- function (nm) unlist(lapply(strsplit(nm, split="_"), function(x) x[3]))
    nm         <- names(eflalo)
    idx.col.w  <- grep('KG', nm) # index columns weight
    idx.col.v  <- grep('EURO', nm)  # index columns value
    v          <- apply(eflalo [,idx.col.v], 2, sum, na.rm=TRUE)
    sp.to.merge.v <- names(v)[v < threshold.value] # e.g. threshold in euros
    sp.to.merge.w <- paste('LE_KG_', get.sp(sp.to.merge.v ), sep='')
    sp.to.merge.w <- sp.to.merge.w[sp.to.merge.w != paste("LE_KG_", code, sep='')]
    sp.to.merge.v <- sp.to.merge.v[sp.to.merge.v != paste("LE_EURO_", code, sep='')]
    OTH.w    <- apply(eflalo[,sp.to.merge.w], 1, sum, na.rm=TRUE)
    OTH.v    <- apply(eflalo[,sp.to.merge.v], 1, sum, na.rm=TRUE)
    eflalo  <- eflalo [, !nm %in% c(sp.to.merge.w, sp.to.merge.v)]
    eflalo[,paste("LE_KG_", code, sep='')] <- replace(eflalo[,paste("LE_KG_", code, sep='')], 
                                                  is.na(eflalo[,paste("LE_KG_", code, sep='')]), 0)
    eflalo[,paste("LE_KG_", code, sep='')]   <- eflalo[,paste("LE_KG_", code, sep='')] + OTH.w
    eflalo[,paste("LE_EURO_", code, sep='')] <- replace(eflalo[,paste("LE_EURO_", code, sep='')], 
                                                  is.na(eflalo[,paste("LE_EURO_", code, sep='')]), 0)
    eflalo[,paste("LE_EURO_", code, sep='')] <- eflalo[,paste("LE_EURO_", code, sep='')] + OTH.v
  return(eflalo)
  }
              