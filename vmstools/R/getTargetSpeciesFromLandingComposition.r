

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#   Francois Bastardie (DTU-Aqua) - Feb 2010 - delivered to the Brest meeting
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# output-based method
# i.e. metier identification using PCA and kmeans clustering on trip- based landing profiles.
# assume that the catch profile perfectly reflects the fisher intention

getTargetSpeciesFromLandingComposition <- function(eflalo, general){

eflalo$tripnum  <- factor(paste(eflalo$VE_REF,".", eflalo$FT_REF, sep='') )
eflalo$logevent <- factor(paste(eflalo$VE_REF, ".", eflalo$FT_REF,".", as.character(eflalo$LE_REC),".",
                           as.character(eflalo$LE_CDAT), ".", as.character(eflalo$LE_GEAR),".",  eflalo$LE_SEQNUM, sep=''))
levels(eflalo$logevent) <- 1: length(levels(eflalo$logevent))


# matter of taste...
if(general$method.pca=="trip")      eflalo$target.trip.pca <- eflalo$tripnum # init
if(general$method.pca=="logevent")  eflalo$target.trip.pca <- eflalo$logevent # init


for (a.gear in levels(eflalo$LE_GEAR)) { # DO THE PCA & CLUSTERING PER FISHERY/GEAR

  cat(paste(a.gear,"\n"))

  idx.in.eflalo <- which(eflalo$LE_GEAR %in% a.gear)
  xx <- eflalo[idx.in.eflalo,] # INPUT
  nm      <- names(eflalo)
  idx.col.w <- grep('KG', nm) # index columns weight
  idx.col.v <- grep('EURO', nm) # index columns value
  dd <- xx [,idx.col.w]
  dd <- as.matrix(dd)
  rownames(dd) <- eflalo[idx.in.eflalo,]$logevent


  if(general$method.pca=="trip"){
  ## AGGREGATE EFLALO TO GET CATCHES PER TRIP INSTEAD OF PER LOG EVENT
  idx.in.eflalo <- which(eflalo$LE_GEAR %in% a.gear)
  xx <- eflalo[idx.in.eflalo,] # INPUT
  nm      <- names(eflalo)
  idx.col.w <- grep('KG', nm) # index columns weight
  idx.col.v <- grep('EURO', nm) # index columns value
  idx.col.tripnum <- which(colnames(eflalo)=='tripnum' )
  dd <- xx [,c(idx.col.v,idx.col.tripnum)]
  dd$tripnum <- factor(dd$tripnum)
  res <- matrix(0, ncol=(ncol(dd)-1), nrow=length(levels(dd$tripnum)))
  # as eflalo is in 'wide' format, the aggregation is a litte bit tricky!
  for(sp in 1: (ncol(dd)-1)){
     res[,sp] <- tapply(dd [,sp], dd$tripnum, sum, na.rm=TRUE)
  }
  colnames(res) <- names(dd[,-ncol(dd)])
  rownames(res) <- as.character(levels(dd$tripnum))
  dd <- res
  }


  ## STEP1- ARBITRARY THRESHOLD i.e. aggregate not well-represented species into a "other" category
  # find the species which are e.g. 100% landings for at least one trip
  zz  <- sweep(dd, c(1), apply(dd, c(1), sum,na.rm=TRUE) ,FUN=quote('/'))
  is.species.targetted.by.some.trips <-  apply(zz,2,max, na.rm=TRUE)
  species.to.keep  <- names(is.species.targetted.by.some.trips[is.species.targetted.by.some.trips>=general$threshold1])
  species.to.merge <- names(is.species.targetted.by.some.trips[is.species.targetted.by.some.trips<general$threshold1])
  dd <- dd[,species.to.keep]
  # collate in a 'other' category
  zz[,'LE_EURO_MZZ']   <- zz[,'LE_EURO_MZZ'] + apply(zz[,species.to.merge],1,sum,na.rm=TRUE)
  sum.kg.per.sp2  <- apply(dd[,species.to.keep],2,sum,na.rm=T)
  zz<- zz [, c(names(sum.kg.per.sp2[order(sum.kg.per.sp2, decreasing=TRUE)]), 'LE_EURO_MZZ') ]
  # remove this 'other' category before doing the pca
  zz <- zz[, !colnames(zz)%in%c('LE_EURO_MZZ')]

  ## STEP2- PCA on SPECIES Landings value per trip or logevent
  library(ade4)
  dataForPCA                 <- replace(zz,is.na(zz) | zz <=0 | zz > 1,0) # debug pca
  res.pca                    <- dudi.pca(data.frame(dataForPCA), center=TRUE, scale=F, scannf=FALSE, nf=5)


  ## STEP3- CLUSTERING using k-means (alternative to 'hclust')
  # the pble with k means is that we have to choose an a priori nb of clusters
  # but if some irrelevant clusters exist, they will be automatically removed when redundant 'main.sp.per.cluster'
  # also kmeans is stochastic, then do iterations...
  for(it in 1: general$nbiter){
     cat(paste("kmeans iter",it,"\n",sep=''))
     NBCLUSTERS <- 20 # choose a nb of cluster large enought
     k     <- kmeans(res.pca$li, NBCLUSTERS, iter.max=10)
    if(FALSE){
      plot(res.pca$li[,"Axis1"], res.pca$li[,"Axis3"], pch="+")
      points(res.pca$li[,"Axis1"],res.pca$li[,"Axis3"], col=k$cluster, pch=3, lwd=3)
    }
    nm                   <- colnames(dataForPCA)
    nm2                  <- paste("kcluster",it,sep='')
    dataForPCA           <-  cbind.data.frame(dataForPCA, k$cluster)
    colnames(dataForPCA) <- c(nm, nm2)
    idx.v                <- grep('EURO', nm) # index columns value
    lst.mean.individuals <-  lapply(split(dataForPCA[,c(nm[idx.v], nm2)], dataForPCA[,nm2]), function(x) apply(x, 2, mean, na.rm=TRUE)*100) # mean individuals per cluster
    mat                  <- matrix(unlist(lst.mean.individuals), nrow=length(unique(dataForPCA[,nm2])), byrow=TRUE)
    mat                  <- mat[,-ncol(mat)]
    colnames(mat)        <- nm[idx.v]
    get.sp               <- function (nm) unlist(lapply(strsplit(nm, split="_"), function(x) x[3]))
    colnames(mat)        <- get.sp (colnames(mat))
    rownames(mat)        <- names(lst.mean.individuals)
    main.sp.per.cluster  <- apply(mat,1,function(x) {  x<- x[order(x, decreasing=TRUE)] ; x[x>general$threshold2]} )  # percent for main species (e.g. >15%) in each cluster
    main.sp.per.cluster  <- lapply(main.sp.per.cluster, function(i) {  i[order(names(i))]} )  # order in alphabetic order to merge e.g. COD-NEP with NEP-COD groups
    main.sp.per.cluster  <- unlist(lapply(main.sp.per.cluster, function(i) paste(names(i),collapse="-")))
    main.sp.per.cluster[main.sp.per.cluster==""] <- "multisp" # in case of clusters with all species e.g. < 15 % of total ladings
    dataForPCA[,nm2]     <- factor(dataForPCA[,nm2])
    levels(dataForPCA[,nm2]) <- main.sp.per.cluster # mapping
    levels(dataForPCA[,nm2])
    table(dataForPCA[,nm2]) # check nb of occurence
   } # end it

 # get the most frequent kmeans among iterations
 nm                  <- colnames(dataForPCA)
 idx.kclust          <- grep('kcluster', nm) # index columns
 dataForPCA$kcluster <- apply(dataForPCA[, idx.kclust], 1, function(x) {cc<- table(x); names(cc[order(cc, decreasing=TRUE)])[1]} )


 # and then: mapping
 tmp <- levels(eflalo$target.trip.pca)[levels(eflalo$target.trip.pca) %in% rownames(dataForPCA)]
 levels(eflalo$target.trip.pca)[levels(eflalo$target.trip.pca) %in% rownames(dataForPCA)] <-
             as.character(dataForPCA[tmp,"kcluster"]) # OUTPUT  : MAPPING
 } # end for a.gear


eflalo$target.trip.pca <- factor(eflalo$target.trip.pca)

cat(paste("target.trip.pca...OK\n",sep=""))

return(eflalo$target.trip.pca)
}

##---------
## CALL
#  eflalo2 <- read.table(file.path("W:","DFAD","Data","eflalo",
#                         paste("eflalo2_",a.year,".csv",sep='')),sep=",", header=TRUE)
# get target and add
#  eflalo2$target.trip.pca <-  getTargetSpeciesFromLandingComposition (eflalo2, general=list(method.pca='trip',
#                                         threshold1=1, threshold2=15, nbiter=7))

