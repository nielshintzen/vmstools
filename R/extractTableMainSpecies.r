################################################################
#   STEP 1 OF THE MULTIVARIATE CLASSIFICATION :                #
#          SELECT MAIN SPECIES OF THE DATASET                  #
#          AND                                                 #
#          EXTRACT A NEW DATASET WITH ONLY THE MAIN SPECIES    #
################################################################


extractTableMainSpecies = function(dat,namesMainSpeciesHAC,paramTotal=95,paramLogevent=100){

    print("######## STEP 1 COMBINATION OF MAIN SPECIES FROM THE THREE EXPLORATORY METHODS ########")
    t1 <- Sys.time()

    # TOTALE
    
    p=ncol(dat)   # Number of species +1
    print("Calculating proportions...")
    propdat=transformation_proportion(dat[,2:p])
    namesSpecies=colnames(propdat)

    # Total quantity caught species by species
    sumcol=rep(as.numeric(NA),p-1) #numeric()
    for(i in 2:p) sumcol[i-1]=sum(dat[,i])
    names(sumcol)=namesSpecies

    # Percent of each species in the total catch
    propSp=sumcol/sum(sumcol,na.rm=TRUE)*100
    # Columns number of each species by decreasing order of capture
    numSp=order(propSp,decreasing=TRUE)
    # Percent of each species in the total catch by decreasing order
    propSp=cumsum(propSp[order(propSp,decreasing=TRUE)])

    Store(objects())
    gc(reset=TRUE)

    # We are taking all species until having at least param1% of total catch
    if (is.null(paramTotal) | !is.numeric(paramTotal)) stop("param1 must be numeric between 0 and 100")
    threshold=paramTotal
    pourcent=which(propSp<=threshold)
    # We are taking the name of selected species
    selSpecies=numSp[1:(length(pourcent)+1)]
    namesSelSpeciesTotal=namesSpecies[selSpecies]



    # LOGEVENT
    
    if (is.null(paramLogevent) | !is.numeric(paramLogevent)) stop("paramLogevent must be numeric between 0 and 100")

    threshold=paramLogevent
    # Selection of species making up over param1% of logevent's captures
    namesSelSpeciesLogevent=character()
    for (i in namesSpecies) {
       if (!is.na(any(propdat[,i]>=threshold)) && any(propdat[,i]>=threshold)) namesSelSpeciesLogevent=c(namesSelSpeciesLogevent,i)
    }

    # Merge with explospecies
    listSpeciesAll <- sort(unique(c(namesMainSpeciesHAC,namesSelSpeciesTotal,namesSelSpeciesLogevent)))
    listSpeciesAll <- listSpeciesAll[!listSpeciesAll=="MZZ"]

    # We are building the table with main species and aggregated other species
    datSpecies=building_tab_pca(propdat,listSpeciesAll)
    rownames(datSpecies)=dat[,1]


    print(" --- end of step 1 ---")
    print(Sys.time()-t1)
    
    return(datSpecies)

}

