# create the VE table to upload in fishframe
# require: data.table package
 mergedTable2FishframeVE <- function (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="DNK")){

    for(what in c('value','weight')){
       load(file.path(general$output.path, paste("all_merged_",what,"_",general$a.year,"_bysn.RData",sep='')))
       nm <- colnames(all.merged)
       if (what=='value') idx.col  <- grep('EURO', nm) # index columns with species
       if (what=='weight') idx.col  <- grep('KG', nm) # index columns with species
       assign(what, apply(all.merged[,idx.col], 1, sum, na.rm=TRUE))
       }
    all.merged$totvalue <- value
    all.merged$totweight <- weight

    an <- function(x) as.numeric(as.character(x))
    all.merged <- all.merged[, !colnames(all.merged) %in% c('fuelcons', 'flag')]
    all.merged <- all.merged[all.merged$SI_STATE==1,] # keep only fishing pings
    nm <- colnames(all.merged)
    all.merged$c_square  <-CSquare(an(all.merged$SI_LONG), an(all.merged$SI_LATI), degrees=0.05)
    all.merged$month <- factor(format(as.POSIXct(all.merged$SI_DATE), "%m"))  # add month
    all.merged$LE_EFF_VMS <- an(all.merged$LE_EFF_VMS) / 24 # convert in hours
    all.merged <- all.merged[,c("LE_EFF_VMS","KW_HOURS","totvalue", "totweight", "LE_MET_level6","c_square","month")]
    all.merged$c_square <- factor(all.merged$c_square)

    # base::aggregate() replaced by fast grouping using the data.table library
    library(data.table)
    DT <- data.table(all.merged)
    qu = quote(list(sum(an(LE_EFF_VMS)),sum(an(KW_HOURS)),sum(an(totvalue)),sum(an(totweight))))
    ff.ve <- DT[,eval(qu), by=list(c_square,month,LE_MET_level6)]
    colnames(ff.ve) <- c('c_square','month','LE_MET_level6','hours','kw_hours', 'totvalue','totweight')

    # additional
    ff.ve$year       <- general$a.year
    ff.ve$country    <- general$a.country
    ff.ve$recordtype <- "VE"


  return(ff.ve)
 }
