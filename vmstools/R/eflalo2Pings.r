#' Merge eflalo to tacsat pings
#' 
#' Internal function of splitAmongPings to merge eflalo landings or values to
#' tacsat pings
#' 
#' 
#' @param eflalo eflalo dataset
#' @param tacsat tacsat dataset
#' @param pings number of pings by variable
#' @param vars variable to merge eflalo to tacsat
#' @param eflaloCol column names of eflalo
#' @param remainTacsat number of tacsat pings that have not been merged yet
#' @return Returns a list of the eflalo dataset, but without the landings and
#' values that have been merged, returns the merged tacsat dataset and returns
#' the number tacsat pings that have not been merged yet
#' @author Niels T. Hintzen
#' @seealso \code{\link{splitAmongPings}}, \code{\link{mergeEflalo2Tacsat}},
#' \code{\link{mergeEflalo2Pings}}
#' @references EU Lot 2 project
#' @export eflalo2Pings
  eflalo2Pings <- function(eflalo,tacsat,pings,vars,eflaloCol,remainTacsat,by=NULL){
        #- Merge landings and values to get unique eflalo set given 'totVars'
      for(iVars in 1:length(vars)){
        eflalo[,vars[iVars]] <- af(ac(eflalo[,vars[iVars]]))
        tacsat[,vars[iVars]] <- af(ac(tacsat[,vars[iVars]]))
      }

      DT                    <- data.table(eflalo)
      eq1                   <- c.listquote(paste("sum(",colnames(eflalo[,kgeur(colnames(eflalo))]),",na.rm=TRUE)",sep=""))
      eq2                   <- c.listquote(vars)

      eflalo            <- data.frame(DT[,eval(eq1),by=eval(eq2)]); colnames(eflalo) <- c(vars,eflaloCol)
      eflalo$ID         <- 1:nrow(eflalo)

        #- Merge eflalo to pings to get number of pings per eflalo record
      byPing             <- merge(eflalo,data.frame(pings),by=vars,all=FALSE)
      byTacsat           <- merge(tacsat,byPing,by=vars,all=FALSE)
      if(is.null(by)==FALSE)
        byTacsat$pings   <- byTacsat$pings / byTacsat[,by]

      try(print(paste("kg in eflalo",round(sum(byPing [,kgeur(colnames(byPing))])))))
      try(print(paste("kg in merged tacsat",round(sum(sweep(byTacsat[,kgeur(colnames(byTacsat))],1,byTacsat$pings,"/"))))))

        #- Bookkeeping which tacsat ID's have been merged and which have not yet been merged
      remainTacsat          <- remainTacsat[which(!remainTacsat %in% byTacsat$ID.x)]

        #- Bookkeeping which eflalo catches have been merged and which have not yet been merged
      idx                   <- sort(unique(byPing$ID))
      try(print(paste("kg removed from eflalo",round(sum(eflalo[idx,kgeur(colnames(eflalo))])))))
      eflalo[idx,kgeur(colnames(eflalo))]   <- 0

    return(list(eflalo=eflalo,tacsat=byTacsat,remainTacsat=remainTacsat))}
