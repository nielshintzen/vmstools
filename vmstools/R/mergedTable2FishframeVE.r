# TO FISHFRAME (author: F. Bastardie)
# create the VE table to upload in fishframe
# required: the data.table package
# optional: the ICES_areas shape file (if not provided as arg then loaded from the vmstools library)
 mergedTable2FishframeVE <- function (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="DNK"),...){
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

    an <- function(x) as.numeric(as.character(x))
    all.merged$SI_LONG <- an( all.merged$SI_LONG)
    all.merged$SI_LATI <- an( all.merged$SI_LATI)
    all.merged <- all.merged[, !colnames(all.merged) %in% c('fuelcons', 'flag')]
    all.merged <- all.merged[all.merged$SI_STATE==1,] # keep only fishing pings
    nm <- colnames(all.merged)
    if(length(lstargs$shape.file)==0) {
    #   data(ICESareas)
    #   all.merged$ICES_area <- ICESarea(long= all.merged$SI_LONG, lat= all.merged$SI_LATI) # DEADLY SLOW!
    all.merged$ICES_area <- ffarea (all.merged[,c('SI_LONG','SI_LATI')])
    }else{
       #all.merged$ICES_area <- ICESarea(long= all.merged$SI_LONG, lat= all.merged$SI_LATI,shape.file=lstargs$shape.file) # DEADLY SLOW!
       all.merged$ICES_area <- ffarea (all.merged[,c('SI_LONG','SI_LATI')])
    }
    all.merged$c_square  <-CSquare(an(all.merged$SI_LONG), an(all.merged$SI_LATI), degrees=0.05)
    all.merged$month <- factor(format(as.POSIXct(all.merged$SI_DATE), "%m"))  # add month
    all.merged$LE_VMS_EFF <- an(all.merged$LE_EFF_VMS) / 24 # convert in hours
    all.merged <- all.merged[,c("LE_EFF_VMS","KW_HOURS","totvalue", "totweight", "LE_MET_level6","ICES_area","c_square","month")]
    all.merged$c_square <- factor(all.merged$c_square)
    all.merged$ICES_area <- factor(all.merged$ICES_area)

 
    # base::aggregate() replaced by fast grouping using the data.table library
    library(data.table)
    DT <- data.table(all.merged)
    qu = quote(list(sum(an(LE_EFF_VMS)),sum(an(KW_HOURS)),sum(an(totvalue)),sum(an(totweight))))
    ff.ve <- DT[,eval(qu), by=list(c_square,ICES_area, month,LE_MET_level6)]
    colnames(ff.ve) <- c('c_square','ICES_area', 'month','LE_MET_level6','hours','kw_hours', 'totvalue','totweight')

    # additional
    ff.ve$year        <- general$a.year
    ff.ve$country     <- general$a.country
    ff.ve$nationalFAC <- NA
    ff.ve$recordtype  <- "VE"
    ff.ve$quarter <- ff.ve$month  # init
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
