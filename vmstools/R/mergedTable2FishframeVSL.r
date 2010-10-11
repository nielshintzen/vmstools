 # TO FISHFRAME
 # create the VSL table to upload in fishframe
 # require: data.table package
 mergedTable2FishframeVSL <- function (general=list(output.path=file.path("C:","output"),
                                          a.year=2009, a.country="DNK")){

      an <- function (x) as.numeric(as.character(x))

      # reshape in 'long' format
         # 1. load
         what <- "weight"
         load(file.path(general$output.path, paste("all_merged_",what,"_",general$a.year,".RData",sep='')))
         nm <- colnames(all.merged)
         idx.col  <- grep('KG', nm) # index columns with species
         all.merged$c_square  <- factor(CSquare(an(all.merged$SI_LONG), an(all.merged$SI_LATI), degrees=0.05))
         all.merged$month <- factor(format(as.POSIXct(all.merged$SI_DATE), "%m"))  # add month
         nm1 <- colnames(all.merged)
         idx.c <- which (nm1 %in% c('VE_REF', 'FT_REF',"LE_MET_level6","c_square","month"))
         xx1 <-  all.merged [, c(idx.c,idx.col)]
         colnames(xx1) <- c('VE_REF', 'FT_REF',"LE_MET_level6","c_square","month", paste( "sp", 1:length(idx.col),sep='') )

         what <- "value"
         load(file.path(general$output.path, paste("all_merged_",what,"_",general$a.year,"_bysn.RData",sep='')))
         nm <- colnames(all.merged)
         idx.col  <- grep('EURO', nm) # index columns with species
         all.merged$c_square  <- factor(CSquare(an(all.merged$SI_LONG), an(all.merged$SI_LATI), degrees=0.05))
         all.merged$month <- factor(format(as.POSIXct(all.merged$SI_DATE), "%m"))  # add month
         nm2 <- colnames(all.merged)
         idx.c <- which (nm2 %in% c('VE_REF', 'FT_REF',"LE_MET_level6","c_square","month"))
         xx2 <-  all.merged [, c(idx.c,idx.col)]
         colnames(xx2) <- c('VE_REF', 'FT_REF',"LE_MET_level6","c_square","month", paste( "sp", 1:length(idx.col),sep='') )

         # 2. reshape
         # (tricky because sub block by sub-block because of 'out of memory')
         res <- NULL
         chunk  <- c( seq(1, nrow(xx1), by=100000), nrow(xx1))
         for(i in 1: (length(chunk)-1)){
            rm(vsl.ff1,vsl.ff2,vsl.ff) ; gc(reset=TRUE)
            cat(paste("lines",chunk[i],"to",chunk[i+1] ,"\n"))
            vsl.ff1 <- reshape( xx1[ chunk[i]:chunk[i+1], ] , direction="long", varying=6:(5+length(idx.col)), sep="") # 'long' format
            colnames(vsl.ff1) <- c('VE_REF', 'FT_REF',"LE_MET_level6","c_square","month", "species", "weight","id")
            vsl.ff1$species <- factor (vsl.ff1$species)
            get.sp <- function (nm) unlist(lapply(strsplit(nm, split="_"), function(x) x[3]))
            levels(vsl.ff1$species) <- get.sp(nm1[idx.col])

            vsl.ff2 <- reshape( xx2[ chunk[i]:chunk[i+1], ] , direction="long", varying=6:(5+length(idx.col)), sep="") # 'long' format
            colnames(vsl.ff2) <- c('VE_REF', 'FT_REF',"LE_MET_level6","c_square","month", "species", "value","id")
            vsl.ff2$species <- factor (vsl.ff2$species)
            nm <- colnames(xx2)
            get.sp <- function (nm) unlist(lapply(strsplit(nm, split="_"), function(x) x[3]))
            levels(vsl.ff2$species) <- get.sp(nm2[idx.col])

            # 3. cbind
            vsl.ff <- cbind.data.frame(vsl.ff1, vsl.ff2$value)

            # 4. clean up
            vsl.ff <- vsl.ff[!is.na(vsl.ff$weight),]
            vsl.ff <- vsl.ff[, !colnames(vsl.ff) %in% c('id')]
            colnames(vsl.ff) <-  c('VE_REF', 'FT_REF',"LE_MET_level6","c_square","month", "species", "weight", "value")

           # 5. aggregate with fast grouping
           library(data.table)
           DT <- data.table(vsl.ff)
           qu = quote(list(sum(an(weight)),sum(an(value))))
           vsl.ff <- DT[,eval(qu), by=list(species,c_square,month,LE_MET_level6)]
           colnames(vsl.ff ) <- c('species','c_square','month','LE_MET_level6','weight','value')

           # 6. bind all
           res <- rbind.data.frame(res, vsl.ff)
           }


   # save
  vsl.ff <- res
  vsl.ff$LE_MET_level6  <- as.character(vsl.ff$LE_MET_level6)
  write.csv2(vsl.ff, file=file.path(general$output.path,
         paste("ff_vsl_", general$a.year, ".csv", sep='')),  row.names = FALSE)


 return(vsl.ff)
 }

      # call
# mergedTable2FishframeVSL (general=list(output.path=file.path("H:","DIFRES","VMSanalysis","results_merged","DKWaters"),
#                                          a.year=2009, a.country="DNK"))
