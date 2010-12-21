

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
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

