 # TO FISHFRAME (author: F. Bastardie)
 # create the VSL table to upload in fishframe
 # require: the 'data.table' and 'doBy' packages
 # optional: the ICES_areas shape file (if not provided as arg then loaded from the vmstools library)
 mergedTable2FishframeVSL <- function (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="DNK"),...){
      lstargs <- list(...)
      
      an <- function (x) as.numeric(as.character(x))
      
      # reshape in 'long' format
         # 1. load
         what <- "weight"
         load(file.path(general$output.path, paste("all_merged_",what,"_",general$a.year,".RData",sep='')))
         nm <- colnames(all.merged)
         idx.col  <- grep('KG', nm) # index columns with species
         all.merged$SI_LONG <- an( all.merged$SI_LONG)
         all.merged$SI_LATI <- an( all.merged$SI_LATI)
         all.merged <- all.merged[all.merged$SI_STATE==1,] # keep only fishing pings
         if(length(lstargs$shape.file)==0) {
         #   data(ICESareas)
         #   all.merged$ICES_area <- ICESarea(long= all.merged$SI_LONG, lat= all.merged$SI_LATI) # DEADLY SLOW!
         all.merged$ICES_area <- get.ICESarea (all.merged[,c('SI_LONG','SI_LATI')])
         }else{
         #all.merged$ICES_area <- ICESarea(long= all.merged$SI_LONG, lat= all.merged$SI_LATI,shape.file=lstargs$shape.file) # DEADLY SLOW!
         all.merged$ICES_area <- get.ICESarea (all.merged[,c('SI_LONG','SI_LATI')])
         }
         all.merged$c_square  <- factor(CSquare(an(all.merged$SI_LONG), an(all.merged$SI_LATI), degrees=0.05))
         all.merged$month <- factor(format(as.POSIXct(all.merged$SI_DATE), "%m"))  # add month
         nm1 <- colnames(all.merged)
         idx.c <- which (nm1 %in% c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area","c_square","month"))
         xx1 <-  all.merged [, c(idx.c,idx.col)]
         colnames(xx1) <- c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area","c_square","month", paste( "sp", 1:length(idx.col),sep='') )
         
         what <- "value"
         load(file.path(general$output.path, paste("all_merged_",what,"_",general$a.year,".RData",sep='')))
         nm <- colnames(all.merged)
         idx.col  <- grep('EURO', nm) # index columns with species
         all.merged$SI_LONG <- an( all.merged$SI_LONG)
         all.merged$SI_LATI <- an( all.merged$SI_LATI)
         all.merged <- all.merged[all.merged$SI_STATE==1,] # keep only fishing pings
         if(length(lstargs$shape.file)==0) {
         #   data(ICESareas)
         #   all.merged$ICES_area <- ICESarea(long= all.merged$SI_LONG, lat= all.merged$SI_LATI) # DEADLY SLOW!
         all.merged$ICES_area <- get.ICESarea (all.merged[,c('SI_LONG','SI_LATI')])
         }else{
         #all.merged$ICES_area <- ICESarea(long= all.merged$SI_LONG, lat= all.merged$SI_LATI,shape.file=lstargs$shape.file) # DEADLY SLOW!
         all.merged$ICES_area <- get.ICESarea (all.merged[,c('SI_LONG','SI_LATI')])
         }
        all.merged$c_square  <- factor(CSquare(an(all.merged$SI_LONG), an(all.merged$SI_LATI), degrees=0.05))
         all.merged$month <- factor(format(as.POSIXct(all.merged$SI_DATE), "%m"))  # add month
         nm2 <- colnames(all.merged)
         idx.c <- which (nm2 %in% c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area","c_square","month"))
         xx2 <-  all.merged [, c(idx.c,idx.col)]
         colnames(xx2) <- c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area","c_square","month", paste( "sp", 1:length(idx.col),sep='') )
       
         # 2. order before splitting in sub-blocks because aggregate() afterwards
         library(doBy)
         xx1 <- orderBy(~c_square+LE_MET_level6+month, data=xx1)
         xx2 <- orderBy(~c_square+LE_MET_level6+month, data=xx2)
         
         # 3. reshape => 'wide' to 'long' format
         # (tricky because sub-block by sub-block because of potential 'out of memory')
         res <- NULL
         lev <- as.character(levels(xx1$c_square)) # do not split randomly but consistently with levels
         chunk  <- c( seq(1, length(lev), by=4000), length(lev)) # 4000 by 4000 levels...
         for(i in 1: (length(chunk)-1)){
            rm(vsl.ff1,vsl.ff2,vsl.ff) ; gc(reset=TRUE)
            cat(paste("level c_square",chunk[i],"to",chunk[i+1] ,"\n"))
            vsl.ff1 <- reshape( xx1[xx1$c_square %in% lev[chunk[i]:chunk[i+1]] , ] , 
                                   direction="long", varying=7:(6+length(idx.col)), sep="") # 'long' format
            colnames(vsl.ff1) <- c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area","c_square","month", "species", "weight","id")
            vsl.ff1$species <- factor (vsl.ff1$species)
            get.sp <- function (nm) unlist(lapply(strsplit(nm, split="_"), function(x) x[3]))
            levels(vsl.ff1$species) <- get.sp(nm1[idx.col])

            vsl.ff2 <- reshape( xx2[ xx2$c_square %in% lev[chunk[i]:chunk[i+1]] , ] ,
                                direction="long", varying=7:(6+length(idx.col)), sep="") # 'long' format
            colnames(vsl.ff2) <- c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area","c_square","month", "species", "value","id")
            vsl.ff2$species <- factor (vsl.ff2$species)
            nm <- colnames(xx2)
            get.sp <- function (nm) unlist(lapply(strsplit(nm, split="_"), function(x) x[3]))
            levels(vsl.ff2$species) <- get.sp(nm2[idx.col])
    
            # 4. cbind
            vsl.ff <- cbind.data.frame(vsl.ff1, vsl.ff2$value) 
         
            # 5. clean up
            vsl.ff <- vsl.ff[!is.na(vsl.ff$weight) & vsl.ff$weight!=0,]
            vsl.ff <- vsl.ff[, !colnames(vsl.ff) %in% c('id')]
            colnames(vsl.ff) <-  c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area", "c_square","month", "species", "weight", "value")
     
           # 6. aggregate with fast grouping 
           library(data.table)
           vsl.ff$ICES_area <- factor(vsl.ff$ICES_area)
           DT <- data.table(vsl.ff)
           qu = quote(list(sum(an(weight)),sum(an(value))))
           vsl.ff <- DT[,eval(qu), by=list(species,ICES_area,c_square,month,LE_MET_level6)]
           colnames(vsl.ff ) <- c('species','ICES_area','c_square','month','LE_MET_level6','weight','value')

           # 7. bind all chunks
           res <- rbind.data.frame(res, vsl.ff)
           }
 
  # 8. additional
  res$year        <- general$a.year
  res$country     <- general$a.country
  res$nationalFAC <- NA
  res$recordtype  <- "VSL"
  res$quarter <- res$month  # init
  levels(res$quarter) <- c(1,1,1,2,2,2,3,3,3,4,4,4)


  # 9. order colums
  ff.vsl <- res
  ff.vsl <- as.data.frame(ff.vsl)[, c('recordtype','country','year','quarter', 'month', 
               'ICES_area','c_square', 'nationalFAC', 'LE_MET_level6',
                 'species','weight','value')]
  # 10. save
  write.table(ff.vsl, file=file.path(general$output.path, 
         paste("ff_vsl_", general$a.year, ".txt", sep='')),  dec=".", sep=";", quote=FALSE, row.names=FALSE)

 
 return(ff.vsl)
 }
}
         
 # call
 #tmp <- mergedTable2FishframeVSL (general=list(output.path=file.path("H:","DIFRES","VMSanalysis","results_merged","DKWaters"),
 #                                         a.year=2009, a.country="DNK"), shape.file=file.path("H:","DIFRES","VMSanalysis","background_map_shape_files","ICES_areas"))
 
