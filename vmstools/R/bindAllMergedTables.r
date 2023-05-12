##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##


#' bind all individual merged tables into a big one
#' 
#' Bind all individual merged tables into one table.
#' 
#' Possibly, merge irrelevant species into an other category to reduce the size
#' of the output object. Possibly, keep only the weight, or the value per
#' species. Both by default, but possible out of memory crash.
#' 
#' @param vessels Vessel ID
#' @param a.year Year of the data
#' @param species.to.keep Array of species to keep
#' @param folder Local folder path
#' @param all.in.one.table Logical: if all data should be merged into one table
#' @author Francois Bastardie
#' @examples
#' 
#' \dontrun{
#'   data(eflalo)
#'   data(tacsat)
#'   data(euharbours); euharbours <- harbours
#' 
#'   # format
#'   eflalo <- formatEflalo(eflalo)
#'   tacsat <- formatTacsat(tacsat)
#' 
#'   # order tacsat chronologically with library(doBy)
#'   tacsat <- sortTacsat(tacsat)
#' 
#'   # test each ping if in harbour or not
#'   tacsat$SI_HARB <- NA
#'   euharbours$Description <- euharbours$harbour
#'   tacsat$SI_HARB <- pointInHarbour(lon=anf(tacsat$SI_LONG),
#'                                    lat=anf(tacsat$SI_LATI),
#'                                    harbours=euharbours,
#'                                    rowSize=30, returnNames=TRUE)
#'   inHarb <- tacsat$SI_HARB
#'   inHarb <- replace(inHarb, !is.na(inHarb), 1)
#'   inHarb <- replace(inHarb, is.na(inHarb), 0)
#'   inHarb <- as.numeric(inHarb)
#' 
#'   # assign a trip identifier
#'   tacsat$SI_FT <- 1 # init
#'   idx <- which(inHarb==0)
#'   tacsat[idx,"SI_FT"] <- cumsum(inHarb) [idx] # add a SI_FT index
#' 
#'   # keep 'out of harbour' points only
#'   # (but keep the departure point and the arrival point lying in the harbour)
#'   startTrip <- c(diff(tacsat[,"SI_FT"]), 0)
#'   endTrip   <- c(0, diff(tacsat[,"SI_FT"]))
#'   tacsat[which(startTrip>0),"SI_FT"]  <-  tacsat[which(startTrip>0)+1,"SI_FT"] 
#'   tacsat[which(endTrip<0),"SI_FT"]    <-  tacsat[which(endTrip<0)-1,"SI_FT"] 
#'   tacsat <- tacsat[which(inHarb==0 |  startTrip>0 |  endTrip<0),]
#' 
#' 
#'   # assign a state to each ping (here, useless if detectFishing at TRUE)
#'   tacsat$SI_STATE <- 2 # init (1: fishing; 2: steaming)
#'   # fake speed rule for fishing state
#'   tacsat$SI_STATE [(tacsat$SI_SP>4 & tacsat$SI_SP<8)] <-1
#' 
#' 
#'   # reduce the size of the eflalo data by merging species
#'   # (assuming that the other species is coded MZZ), threshold in euros.
#'   eflalo2 <- poolEflaloSpecies (eflalo, threshold=1e6, code="MZZ")
#' 
#'   # debug if eflalo has not been cleaned earlier
#'   eflalo <- eflalo[!eflalo$VE_REF=="NA" &!is.na(eflalo$VE_REF),]
#'   
#'   # an informed VE_FLT is also required
#'   if(all(is.na(eflalo$VE_FLT))) eflalo$VE_FLT <- "fleet1"
#'   
#'   # possible mis-naming mistakes
#'     if(!match('LE_MET_level6',colnames(eflalo))>0){
#'       eflalo$LE_MET_level6 <- eflalo$LE_MET
#'     }
#' 
#'   # debug
#'   eflalo <- eflalo[eflalo$LE_MET!="No_logbook6",]
#' 
#' 
#'   # TEST FOR A GIVEN SET OF VESSELS
#'   # (if detect.fishing is true then do also detection of fishing activity
#'   # e.g. if speed='segment' the segmentTacsatSpeed() automatic detection of fishing states
#'   # that will overwrite the existing SI_STATE)
#'   mergeEflalo2Pings (eflalo=eflalo, tacsat=tacsat, vessels=c("738", "804"),
#'                      general=list(output.path=file.path("C:","output"),
#'                      visual.check=TRUE, detectFishing=TRUE, speed="segment",
#'                      what.speed="calculated"))
#'   # ...OR APPLY FOR ALL VESSELS IN eflalo
#'   mergeEflalo2Pings (eflalo=eflalo, tacsat=tacsat,
#'                      general=list(output.path=file.path("C:","output"),
#'                      visual.check=TRUE, detectFishing=TRUE, speed="segment",
#'                      what.speed="calculated"))
#'   gc(reset=TRUE)
#' 
#'   # load the merged output table for one vessel
#'   load(file.path("C:","output","merged_804_1800.RData"))
#' 
#'   # check the conservation of landings
#'   sum(tapply(anf(merged$LE_KG_PLE), merged$flag, sum, na.rm=TRUE))
#'   sum(eflalo[eflalo$VE_REF=="804","LE_KG_PLE"], na.rm=TRUE)
#' 
#' 
#'    # ...or bind all vessels (keeping only some given species here)
#'   bindAllMergedTables (vessels=c("738", "804"), a.year = "1800",
#'                       species.to.keep=c("PLE","COD"),
#'                       folder = file.path("C:","output"),
#'                       all.in.one.table=TRUE)
#' 
#'     # ...and load the merged output table for all vessels
#'   load(file.path("C:","output","all_merged__1800.RData"))
#' 
#' }
#' 
#' 
#' @export bindAllMergedTables
  bindAllMergedTables <- 
    function (vessels=character(), a.year='2009',species.to.keep=character(), 
                      folder = file.path(), all.in.one.table=FALSE){

  
      bindAll <- function(vessels, a.year, species.to.keep, what, folder){
                    
        all.merged <- NULL # init
        count <- 0
        for(a.vesselid in as.character(vessels)){
          count <- count+1
          print(count/length(vessels)*100)
          cat(paste("load 'merged' for",a.vesselid,"\n"))
          er <- try(load(file = file.path(folder, paste('merged_',a.vesselid,'_',a.year,'.RData',sep=''))) 
                        , silent=TRUE) # get the 'merged' table for this vessel

 
         if(class(er)!="try-error")   {
           if(length(species.to.keep)!=0){
              get.sp            <- function (nm) unlist(lapply(strsplit(nm, split="_"), function(x) x[3]))
              nm                <- names(merged)
              idx.col.w         <- grep('KG', nm) # index columns weight
              all.sp            <- get.sp(nm[idx.col.w])
              species.to.merge <- all.sp[!all.sp %in% species.to.keep]
              # merge to other sp
              merged$LE_EURO_MZZ <-  replace(merged$LE_EURO_MZZ, is.na(merged$LE_EURO_MZZ), 0)
              merged$LE_EURO_MZZ <-  merged$LE_EURO_MZZ + apply(merged[, paste('LE_EURO_',species.to.merge,sep='')], 1, sum, na.rm=TRUE)
              merged <-  merged[, !colnames(merged) %in% paste('LE_EURO_',species.to.merge,sep='')]
              merged$LE_KG_MZZ <-  replace(merged$LE_KG_MZZ, is.na(merged$LE_KG_MZZ), 0)
              merged$LE_KG_MZZ <-  merged$LE_KG_MZZ + apply(merged[, paste('LE_KG_',species.to.merge,sep='')], 1, sum, na.rm=TRUE)
              merged <-  merged[, !colnames(merged) %in% paste('LE_KG_',species.to.merge,sep='')]
              }
           nm       <- names(merged)
           idx.col.w <- grep('KG', nm) # index columns with species weight
           idx.col.v <- grep('EURO', nm) # index columns with species value
           if(length(what)==0){  idx.col   <- c(idx.col.w, idx.col.v) # get all (but possible out of memory crash)
                  } else{ 
                    if(what=='weight') idx.col   <- idx.col.w
                    if(what=='value')  idx.col   <- idx.col.v
                  }
           # KEEP THE ESSENTIAL
           merged <- merged[,c( "VE_REF", "FT_REF", "VE_FLT", "LE_MET_level6", "LE_GEAR",
                                  "SI_LATI","SI_LONG", "SI_SP", "SI_HE", "SI_STATE", "SI_DATE", "SI_TIME", "SI_HARB",
                                     nm[idx.col], 'LE_EFF_VMS', 'KW_HOURS',
                                         "flag")]
           all.merged <- rbind.data.frame(all.merged, merged)
         } else cat(paste("failure for",a.vesselid,"\n"))
      print(nrow(all.merged))
      }
      # save
      save("all.merged", file = file.path(folder, paste('all_merged_',what,'_',a.year,'.RData',sep='') ))

    return()
    }
    
    # calls
    bindAll (vessels, a.year, species.to.keep, what='weight', folder)
    bindAll (vessels,a.year, species.to.keep, what='value', folder)
    if(all.in.one.table) bindAll (vessels, a.year, species.to.keep, what=character(), folder)



 return()
 }

