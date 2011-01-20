################################################################
#   STEP 1 OF THE MULTIVARIATE CLASSIFICATION :                #
#          SELECT MAIN SPECIES OF THE DATASET                  #
#          AND                                                 #
#          EXTRACT A NEW DATASET WITH ONLY THE MAIN SPECIES    #
################################################################


extractTableMainSpecies = function(dat,NamesMainSpeciesHAC,paramTotal=95,paramLogevent=100){

    print("######## STEP 1 COMBINATION OF MAIN SPECIES FROM THE THREE EXPLORATORY METHODS ########")


    t1 <- Sys.time()

    # clu : first select the appropriate columns
    #dat <- dat[,c("LE_ID",grep(Val,names(dat),value=T))]

    #dat[is.na(dat)] <- 0


    # then simplify the names of columns
    #names(dat)[-1] <- unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))
    #names(dat)

    #list Species Totale
    p=ncol(dat)   # Number of species +1
    #for (i in 2:p) dat[is.na(dat[,i]),i] <-0


    print("calculating proportions...")

    propdat=transformation_proportion(dat[,2:p])
    nameSpecies=colnames(propdat)

    # Total quantity caught species by species
    sumcol=rep(as.numeric(NA),p-1) #numeric()
    for(i in 2:p) sumcol[i-1]=sum(dat[,i])

    names(sumcol) <- names(dat)[-1]

    #sumcol <- apply(dat[,-1],2,sum)

    # Total quantity caught
    # sumtotale=sum(sumcol,na.rm=T)
    # Percent of each species in the total catch
    propesp=sumcol/sum(sumcol,na.rm=T)*100
    # Columns number of each species by decreasing order of capture
    numesp=order(propesp,decreasing=T)
    # Percent of each species in the total catch by decreasing order
    propesp=cumsum(propesp[order(propesp,decreasing=T)])

    Store(objects())#[-which(objects() %in% c('dat','methSpecies','param1','param2','pcaYesNo','methMetier','param3','param4'))])
    gc(reset=TRUE)

    # We are taking all species until having at least param1% of total catch
    if (is.null(paramTotal) | !is.numeric(paramTotal)) stop("param1 must be numeric between 0 and 100")
    seuil=paramTotal
    pourcent=which(propesp<=seuil)
    # We are taking the name of selected species
    espsel=numesp[1:(length(pourcent)+1)]
    nomespselTotal=nameSpecies[espsel]

    #    #eventually removing MZZ
    #    nomespselTotal <- nomespsel[!nomespsel=="MZZ"]
    #    #cat("main species:",sort(nomespsel),"\n")

    #LogEvent
    if (is.null(paramLogevent) | !is.numeric(paramLogevent)) stop("paramLogevent must be numeric between 0 and 100")

    seuil=paramLogevent
    # Selection of species making up over param1% of logevent's captures
#    pourcent <- apply(propdat,1,function(x) which(x>=seuil))
#    nomespselLogevent <- names(propdat)[unique(unlist(pourcent))]
    nomespselLogevent=character()
    for (i in nameSpecies) if (any(propdat[,i]>=seuil)) nomespselLogevent=c(nomespselLogevent,i)

    #eventually removing MZZ
   # nomespselLogevent <- nomespsel[!nomespsel=="MZZ"]

    #merge with explospecies
    ListSpeciesAll <- sort(unique(c(NamesMainSpeciesHAC,nomespselTotal,nomespselLogevent)))
    ListSpeciesAll <- ListSpeciesAll[!ListSpeciesAll=="MZZ"]

    # We are bulding the table with main species and aggregated other species
    datSpecies=building_tab_pca(propdat,ListSpeciesAll)

    #datSpecies <- data.frame(LE_ID=dat$LE_ID,datSpecies)
    rownames(datSpecies)=dat$LE_ID


    print(" --- end of step 1 ---")
    print(Sys.time()-t1)


    return(datSpecies)

}

