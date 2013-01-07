

 pings2Fishframe <- function(general=list(output.path=
             file.path("H:","DIFRES","VMSanalysis","results_merged","DKWaters"),
                                          a.year=2009, a.country="DNK", degree=0.05)){
                                          
  
 # TO FISHFRAME (author: F. Bastardie)
 # create the VE table to upload in fishframe
 # required: the data.table package
 # optional: the  "areas" shape file for ICESarea()(if not provided as arg then loaded from the vmstools library)
 mergedTable2FishframeVE <- function (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="DNK", degree=0.05),...){
    lstargs <- list(...)

    for(what in c('value','weight')){
       load(file.path(general$output.path, paste("all_merged_",what,"_",general$a.year,".RData",sep='')))
       nm <- colnames(all.merged)
       if (what=='value') idx.col  <- grep('EURO', nm) # index columns with species
       if (what=='weight') idx.col  <- grep('KG', nm) # index columns with species
       assign(what, apply(all.merged[,idx.col], 1, sum, na.rm=TRUE))
       }
    all.merged$totvalue <- value
    all.merged$totweight <- weight

    # debug
    all.merged[all.merged$LE_MET_level6=="", "LE_MET_level6"] <-"NA"

    an <- function(x) as.numeric(as.character(x))
    all.merged$SI_LONG <- an( all.merged$SI_LONG)
    all.merged$SI_LATI <- an( all.merged$SI_LATI)
    all.merged <- all.merged[, !colnames(all.merged) %in% c('fuelcons', 'flag')]
    all.merged <- all.merged[all.merged$SI_STATE==1,] # keep only fishing pings
    nm <- colnames(all.merged)
    if(length(lstargs$spatialPolygons)==0) {
       data(ICESareas)
       all.merged$ICES_area <- ICESarea (all.merged[,c('SI_LONG','SI_LATI')], areas=ICESareas)
    }else{
       all.merged$ICES_area <- ICESarea (all.merged[,c('SI_LONG','SI_LATI')], areas=lstargs$spatialPolygons)
    }
    all.merged$c_square   <- CSquare(an(all.merged$SI_LONG), an(all.merged$SI_LATI), degrees=general$degree)
    all.merged$month      <- factor(format(as.POSIXct(all.merged$SI_DATE), "%m"))  # add month
    all.merged$LE_EFF_VMS <- an(all.merged$LE_EFF_VMS) / 24 # convert in hours
    all.merged            <- all.merged[,c("LE_EFF_VMS","KW_HOURS","totvalue", "totweight", "LE_MET_level6","ICES_area","c_square","month")]
    all.merged$c_square   <- factor(all.merged$c_square)
    all.merged$ICES_area  <- factor(all.merged$ICES_area)

 
    # base::aggregate() replaced by fast grouping using the data.table library
    library(data.table)
    DT     <- data.table(all.merged)
    qu     <- quote(list(sum(an(LE_EFF_VMS)),sum(an(KW_HOURS)),sum(an(totvalue)),sum(an(totweight))))
    ff.ve  <- DT[,eval(qu), by=list(c_square,ICES_area, month,LE_MET_level6)]
    colnames(ff.ve) <- c('c_square','ICES_area', 'month','LE_MET_level6','hours','kw_hours', 'totvalue','totweight')

    # additional
    ff.ve$year            <- general$a.year
    ff.ve$country         <- general$a.country
    ff.ve$nationalFAC     <- " "
    ff.ve$recordtype      <- "VE"
    ff.ve$quarter         <- ff.ve$month  # init
    levels(ff.ve$quarter) <- c(1,1,1,2,2,2,3,3,3,4,4,4)


  
  #order colums
  ff.ve <- as.data.frame(ff.ve)[, c('recordtype','country','year','quarter', 'month', 
               'ICES_area','c_square', 'nationalFAC', 'LE_MET_level6',
                 'hours','kw_hours','totweight','totvalue')]
  # save
  write.table(ff.ve, file=file.path(general$output.path, 
         paste("ff_ve_", general$a.year, ".txt", sep='')),  dec=".", sep=";", quote=FALSE, row.names=FALSE)

  return(ff.ve)
 }

 # TO FISHFRAME (author: F. Bastardie)
 # create the VSL table to upload in fishframe
 # require: the 'data.table' and 'doBy' packages
 # optional: the  "areas" shape file for ICESarea()(if not provided as arg then loaded from the vmstools library)
 mergedTable2FishframeVSL <- function (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="DNK", degree=0.05),...){
      lstargs <- list(...)
      
      an <- function (x) as.numeric(as.character(x))
      
      # reshape in 'long' format
         # 1. load
         what <- "weight"
         load(file.path(general$output.path, paste("all_merged_",what,"_",general$a.year,".RData",sep='')))
         all.merged[all.merged$LE_MET_level6=="", "LE_MET_level6"] <-"NA"  # debug
         nm <- colnames(all.merged)
         idx.col  <- grep('KG', nm) # index columns with species
         all.merged$SI_LONG <- an( all.merged$SI_LONG)
         all.merged$SI_LATI <- an( all.merged$SI_LATI)
         all.merged <- all.merged[all.merged$SI_STATE==1,] # keep only fishing pings
         if(length(lstargs$spatialPolygons)==0) {
         all.merged$ICES_area <- ICESarea (all.merged[,c('SI_LONG','SI_LATI')], areas=ICESareas)
         }else{
         data(ICESareas)
         all.merged$ICES_area <- ICESarea (all.merged[,c('SI_LONG','SI_LATI')], areas=lstargs$spatialPolygons)
         }
         all.merged$c_square  <- factor(CSquare(an(all.merged$SI_LONG), an(all.merged$SI_LATI), degrees=general$degree))
         all.merged$month     <- factor(format(as.POSIXct(all.merged$SI_DATE), "%m"))  # add month
         nm1                  <- colnames(all.merged)
         idx.c                <- which (nm1 %in% c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area","c_square","month"))
         xx1                  <-  all.merged [, c(idx.c,idx.col)]
         colnames(xx1)        <- c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area","c_square","month", paste( "sp", 1:length(idx.col),sep='') )
         
         what <- "value"
         load(file.path(general$output.path, paste("all_merged_",what,"_",general$a.year,".RData",sep='')))
         all.merged[all.merged$LE_MET_level6=="", "LE_MET_level6"] <-"NA"  # debug
         nm                 <- colnames(all.merged)
         idx.col            <- grep('EURO', nm) # index columns with species
         all.merged$SI_LONG <- an( all.merged$SI_LONG)
         all.merged$SI_LATI <- an( all.merged$SI_LATI)
         all.merged         <- all.merged[all.merged$SI_STATE==1,] # keep only fishing pings
         if(length(lstargs$spatialPolygons)==0) {
             data(ICESareas)
             all.merged$ICES_area <- ICESarea (all.merged[,c('SI_LONG','SI_LATI')], areas=ICESareas)
         }else{
             all.merged$ICES_area <- ICESarea (all.merged[,c('SI_LONG','SI_LATI')], areas=lstargs$spatialPolygons)
         }
        all.merged$c_square  <- factor(CSquare(an(all.merged$SI_LONG), an(all.merged$SI_LATI), degrees=general$degree))
        all.merged$month     <- factor(format(as.POSIXct(all.merged$SI_DATE), "%m"))  # add month
        nm2                  <- colnames(all.merged)
        idx.c                <- which (nm2 %in% c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area","c_square","month"))
        xx2                  <- all.merged [, c(idx.c,idx.col)]
        colnames(xx2)        <- c('VE_REF', 'FT_REF',"LE_MET_level6","ICES_area","c_square","month", paste( "sp", 1:length(idx.col),sep='') )
       
         # 2. order before splitting in sub-blocks because aggregate() afterwards
         library(doBy)
         xx1 <- orderBy(~c_square+LE_MET_level6+month, data=xx1)
         xx2 <- orderBy(~c_square+LE_MET_level6+month, data=xx2)
         
         # 3. reshape => 'wide' to 'long' format
         # (tricky because sub-block by sub-block because of potential 'out of memory')
         res <- NULL
         lev <- as.character(levels(xx1$c_square)) # do not split randomly but consistently with levels
         chunk  <- c( seq(1, length(lev), by=2000), length(lev)) # 2000 by 2000 levels...
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
    
           # 6. aggregate with fast grouping (caution: > R.2.11.0) 
           library(data.table)
           vsl.ff$ICES_area <- factor(vsl.ff$ICES_area)
           DT      <- data.table(vsl.ff)
           qu      <- quote(list(sum(an(weight)),sum(an(value))))
           vsl.ff  <- DT[,eval(qu), by=list(species,ICES_area,c_square,month,LE_MET_level6)]
           colnames(vsl.ff ) <- c('species','ICES_area','c_square','month','LE_MET_level6','weight','value')

           # 7. bind all chunks
           res <- rbind.data.frame(res, vsl.ff)
           }
 
  # 8. additional
  res$year            <- general$a.year
  res$country         <- general$a.country
  res$nationalFAC     <- " "
  res$recordtype      <- "VSL"
  res$quarter         <- res$month  # init
  levels(res$quarter) <- c(1,1,1,2,2,2,3,3,3,4,4,4)

  # 9. convert species fao code to fishframe latin species names
  data(speciesLatinNames)
  res$species <-  speciesLatinNames$ff_species_latin[match(as.character(res$species),
                               as.character(speciesLatinNames$fao_code))]

   
  # 10. order colums
  ff.vsl <- res
  ff.vsl <- as.data.frame(ff.vsl)[, c('recordtype','country','year','quarter', 'month', 
               'ICES_area','c_square', 'nationalFAC', 'LE_MET_level6',
                 'species','weight','value')]
  # 11. save
  write.table(ff.vsl, file=file.path(general$output.path, 
         paste("ff_vsl_", general$a.year, ".txt", sep='')),  dec=".", sep=";", quote=FALSE, row.names=FALSE)

 
 return(ff.vsl)
 }
         
 
  
  
  # GENERAL CALLS
  ve   <- mergedTable2FishframeVE  (general=general)
  vsl  <- mergedTable2FishframeVSL (general=general)


  # add a fake column to get the same ncol()
  vsl <- cbind(vsl, 0)
  colnames(ve)  <- paste('col', 1:ncol(ve), sep='')
  colnames(vsl) <- paste('col', 1:ncol(vsl), sep='')

  # bind and order
  #(to get a VE line and then VSL lines, VE and then VSL lines, etc.)
  ff <- rbind(ve,vsl)
  library(doBy)
  ff <- orderBy(~col7+col9+col5+col6+col1, data=ff)

  # round the numbers
  ff[, c(11,12,13)] <- ceiling(ff[, c(11,12,13)])

  # remove record with 0 in value because both
  # weight and value are mandatory for uploading in fishframe
  # but this should not have any effect here because 
  # filling the gap should have been performed before on eflalo...
  ff <- ff[ff[, c(12)]!=0,]

  # save
  write.table(ff, file=file.path(general$output.path,
         paste(general$a.country,"_", general$a.year, "VD.csv", sep='')),  dec=".", sep=";",
          quote=FALSE, row.names=FALSE, col.names=FALSE)


  return(ff)
  }


# example calls
# vsl <- mergedTable2FishframeVSL (general=list(output.path=file.path("C:","merging", "EflaloAndTacsat"),
#                                          a.year=2009, a.country="DNK", degree=0.05) )
# ve <- mergedTable2FishframeVE  (general=list(output.path=file.path("C:","merging", "EflaloAndTacsat"),
#                                          a.year=2009, a.country="DNK", degree=0.05) )

# alternatively:
#for (a_year in as.character(2005:2010)) 
#  ff <- pings2Fishframe (general=list(output.path=file.path("C:","merging", "EflaloAndTacsat"),
#                                          a.year=a_year, a.country="DNK", degree=0.01) )

