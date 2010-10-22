# Exploration of Step 1 : SELECTION OF SPECIES

ExploSpeciesSelection=function(dat,analysisName="",Val="EURO"){

dat <- dat[,c("LE_ID",grep(Val,names(dat),value=T))]
dat[is.na(dat)] <- 0  


    # First simplify the names of columns
    names(dat)[-1] <- unlist(lapply(strsplit(names(dat[,-1]),"_"),function(x) x[[3]]))

    # HAC
    print("######## STEP 1 SELECTION OF MAIN SPECIES BY 'HAC' METHOD ########")
    
    toutfait=FALSE
    p=ncol(dat)   # Number of species
    n=nrow(dat)
    # Transform quantities to proportions of total quantity caught by logevent
    print("calculating proportions...") 

    propdat=transformation_proportion(dat[,2:p])
    nameSpecies=names(propdat)
    nameSpecies
    # Transposing data
    table_var=table_variables(propdat)

    # HAC
    print("cluster...")
    cah_var=hcluster(table_var, method="euclidean", link="ward")

    # Select the number of clusters by scree-test
    inerties.vector=cah_var$height[order(cah_var$height,decreasing=T)]
    nb.finalclusters=which(scree(inerties.vector)$epsilon<0)[1]

    # Dendogram cutting at the selected level
    cah_cluster_var=cutree(cah_var,k=nb.finalclusters)

    # Selection of main species
    nomespsel=select_species(dat[,2:p],cah_cluster_var)
    espdata=names(dat[,2:p])
    cat("main species:",nomespsel,"\n")


print("selection of possible residual species...")        
    
    while(toutfait==FALSE){
      ind_princ=which(is.element(espdata,nomespsel))
      ind_autres=setdiff(espdata,nomespsel)
      princ=table_var[ind_princ,]
      resi=table_var[ind_autres,]

      # HAC on residual species group
      resicah_var=hcluster(resi, method="euclidean", link="ward")
      # Select the number of clusters by scree-test
      resiinerties.vector=resicah_var$height[order(resicah_var$height,decreasing=T)]
      resinbclusters=which(scree(resiinerties.vector)$epsilon<0)[1]
      # Cut the dendogram at the selected level
      resicah_cluster_var=cutree(resicah_var,k=resinbclusters)

      # Number and names of species by cluster, mean of capture by species for each cluster
      long=numeric()
      somclus=numeric()
      moyclus=numeric()
      for(i in 1:resinbclusters){
        long[i]=length(which(resicah_cluster_var==i))
        nomclus=names(which(resicah_cluster_var==i))
        somclus[i]=sum(resi[nomclus,])
        moyclus[i]=somclus[i]/long[i]
      }
      # we are taking the cluster which represents the greatest part of capture
      indclusmax=which.max(moyclus)
      names(which(resicah_cluster_var==indclusmax))
      # if the cluster owns only one species, we are adding the species to the main species and restart the CAH with residuals species
      if(length(which(resicah_cluster_var==indclusmax))==1){
        nomespsel=c(nomespsel,names(which(resicah_cluster_var==indclusmax)))
        print(paste("adding species", names(which(resicah_cluster_var==indclusmax)))) 

        toutfait=FALSE
      }else{   # else we are stopping here and keeping the main species of the start
        toutfait=TRUE
      }
    }

    # Return the dataset retaining only the main species
    #datSpecies=building_tab_pca(propdat,nomespsel)
    datSpeciesWithoutProp=building_tab_pca(dat[,2:p],nomespsel)
    nbMainSpeciesHAC=length(nomespsel)
    NamesMainSpeciesHAC <- nomespsel

    pourcentCatchMainSpeciesHAC=apply(datSpeciesWithoutProp,1,sum)/apply(dat[,2:p],1,sum)*100
    #Delete logevents without catch
    logWithoutCatch=which(is.na(pourcentCatchMainSpeciesHAC))
    logWithCatch=setdiff(1:n,logWithoutCatch)

    pourcentCatchMainSpeciesHAC=apply(datSpeciesWithoutProp[logWithCatch,],1,sum)/apply(dat[logWithCatch,2:p],1,sum)*100
    medianPourcentCatchMainSpeciesHAC=median(pourcentCatchMainSpeciesHAC)
    mean(pourcentCatchMainSpeciesHAC)

    
    # TOTALE
    
    print("######## STEP 2 SELECTION OF MAIN SPECIES BY 'TOTALE' METHOD ########")
    
    # Total quantity caught by species
    sumcol=numeric()
#    for(i in 2:p){
#      sumcol[i]=sum(dat[,i])
#    }
    
    sumcol <- apply(dat[,-1],2,sum)

    # Total quantity caught
    sumtotale=sum(sumcol,na.rm=T)
    # Percent of each species in the total catch
    propesp=sumcol/sumtotale*100
    # by decreasing order
    propespdec=propesp[order(propesp,decreasing=T)]
    # Columns number of each species by decreasing order of capture
    numesp=order(propesp,decreasing=T)
    # Cumulative percent of catch
    propespcum=cumsum(propespdec)
    # We are taking all species until having at least seuil% of total catch

    nbMainSpeciesTotale=numeric()
    medianPourcentCatchMainSpeciesTotale=numeric()

    for(seuil in seq(5,100,5)){
      cat("seuil:",seuil,"\n")
      pourcent=which(propespcum<=seuil)
      espsel=numesp[1:(length(pourcent)+1)]
      # We are taking the name of selected species
      nomespsel=nameSpecies[espsel]
      # We are bulding the table with main species and aggregated other species
      #datSpecies=building_tab_pca(propdat,nomespsel)
      datSpeciesWithoutProp=building_tab_pca(dat[,2:p],nomespsel)
      nbMainSpeciesTotale[seuil/5]=length(nomespsel)

      pourcentCatchMainSpeciesTotale=apply(datSpeciesWithoutProp,1,sum)/apply(dat[,2:p],1,sum)*100
      #Delete logevents without catch
      logWithoutCatch=which(is.na(pourcentCatchMainSpeciesTotale))
      logWithCatch=setdiff(1:n,logWithoutCatch)
      if(length(nomespsel)==1){
        vectorNul=rep(0,n)
        datSpeciesWithoutProp=cbind(datSpeciesWithoutProp,vectorNul)
      }
      pourcentCatchMainSpeciesTotale=apply(datSpeciesWithoutProp[logWithCatch,],1,sum)/apply(dat[logWithCatch,2:p],1,sum)*100
      medianPourcentCatchMainSpeciesTotale[seuil/5]=median(pourcentCatchMainSpeciesTotale)
      mean(pourcentCatchMainSpeciesTotale)
    }
    nbMainSpeciesTotale=c(0,nbMainSpeciesTotale)
    medianPourcentCatchMainSpeciesTotale=c(0,medianPourcentCatchMainSpeciesTotale)
    NamesMainSpeciesTotale <- nomespsel[1:nbMainSpeciesHAC]


    # LOGEVENT
    
    print("######## STEP 3 SELECTION OF MAIN SPECIES BY 'LOGEVENT' METHOD ########")
    
    nbMainSpeciesLogevent=numeric()
    medianPourcentCatchMainSpeciesLogevent=numeric()

    #AllNamesSpeciesLogevent <- vector()
    for(seuil in seq(5,100,5)){
    cat("seuil:",seuil,"\n")
      # Selection of species making up over seuil% of logevent's captures
      pourcent <- apply(propdat,1,function(x) which(x>=seuil))
      nomespsel <- names(propdat)[unique(unlist(pourcent))]
      #AllNamesSpeciesLogevent <- c(AllNamesSpeciesLogevent,nomespsel)

      # We are bulding the table with main species and aggregated other species
      #datSpecies=building_tab_pca(propdat,nomespsel)
      datSpeciesWithoutProp=building_tab_pca(dat[,2:p],nomespsel)
      nbMainSpeciesLogevent[seuil/5]=length(nomespsel)
      
      pourcentCatchMainSpeciesLogevent=apply(datSpeciesWithoutProp,1,sum)/apply(dat[,2:p],1,sum)*100
      #Delete logevents without catch
      logWithoutCatch=which(is.na(pourcentCatchMainSpeciesLogevent))
      logWithCatch=setdiff(1:n,logWithoutCatch)
      if(length(nomespsel)==1){
        vectorNul=rep(0,n)
        datSpeciesWithoutProp=cbind(datSpeciesWithoutProp,vectorNul)
      }
      pourcentCatchMainSpeciesLogevent=apply(datSpeciesWithoutProp[logWithCatch,],1,sum)/apply(dat[logWithCatch,2:p],1,sum)*100
      medianPourcentCatchMainSpeciesLogevent[seuil/5]=median(pourcentCatchMainSpeciesLogevent)
    }
    nbMainSpeciesLogevent=c(p-1,nbMainSpeciesLogevent)
    medianPourcentCatchMainSpeciesLogevent=c(100,medianPourcentCatchMainSpeciesLogevent)
    #NamesMainSpeciesLogevent <- names(sort(table(AllNamesSpeciesLogevent),decreasing=T)[1:nbMainSpeciesHAC])
    NamesMainSpeciesLogevent <- nomespsel
    


    # GRAPHICS
    # Number of main species
    png(paste(analysisName,"Number of main species.png",sep="_"), width = 1200, height = 800)
    plot(seq(0,100,5),nbMainSpeciesTotale,type='l',col="blue",lwd=2, main="Number of main species selected depending of the threshold", xlab="Threshold (%)",ylab="Number of species")
    lines(seq(0,100,5),nbMainSpeciesLogevent,col="green",lwd=2)
    abline(nbMainSpeciesHAC,0, col="red",lwd=2)
    mtext(paste(p-1," Species"),col='darkblue')
    legend(70, p*0.9, c( "HAC", "Totale", "Logevent"),lwd=2,col=c("red", "blue", "green"))
    dev.off()

    # Median percentage of catch represented by main species by logevent
    png(paste(analysisName,"Median percentage of catch represented by main species by logevent.png",sep="_"), width = 1200, height = 800)
    plot(seq(0,100,5),medianPourcentCatchMainSpeciesTotale,type='l',col="blue",lwd=2, main="Median percentage of catch represented by main species by logevent depending of the threshold", xlab="Threshold (%)",ylab="Median percentage of catch represented by main species by logevent")
    lines(seq(0,100,5),medianPourcentCatchMainSpeciesLogevent,col="green",lwd=2)
    abline(medianPourcentCatchMainSpeciesHAC,0, col="red",lwd=2)
    mtext(paste(p-1," Species"),col='darkblue')
    legend(70, 40, c("HAC", "Totale", "Logevent"),lwd=2,col=c("red", "blue", "green"))
    dev.off()



    return(list(nbMainSpeciesHAC=nbMainSpeciesHAC, nbMainSpeciesTotale=nbMainSpeciesTotale, nbMainSpeciesLogevent=nbMainSpeciesLogevent,
    NamesMainSpeciesHAC=sort(NamesMainSpeciesHAC), NamesMainSpeciesTotale=sort(NamesMainSpeciesTotale), NamesMainSpeciesLogevent=sort(NamesMainSpeciesLogevent),
    medianPourcentCatchMainSpeciesHAC=medianPourcentCatchMainSpeciesHAC, medianPourcentCatchMainSpeciesTotale=medianPourcentCatchMainSpeciesTotale,
    medianPourcentCatchMainSpeciesLogevent=medianPourcentCatchMainSpeciesLogevent))
    
}
