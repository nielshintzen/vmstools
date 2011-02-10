################################################################################
# EXPLORE THE SELECTED SPECIES DEPENDING ON THE METHOD (HAC, TOTALE, LOGEVENT) #
# AND THE THRESHOLD CHOSEN                                                     #
################################################################################

selectMainSpecies=function(dat,analysisName="",RunHAC=TRUE,DiagFlag=FALSE){

    p=ncol(dat)   # Number of species
    n=nrow(dat)
    
    # Delete logevents without catch  - THIS SHOULD BE DONE OUTSIDE 
    #logWithoutCatch=which(apply(dat[,2:p],1,sum)==0)
    #logWithCatch=setdiff(1:n,logWithoutCatch)
    #dat=dat[logWithCatch,]
    #n=nrow(dat)
      
    # Transform quantities to proportions of total quantity caught by logevent
    print("calculating proportions...") 

    propdat=transformation_proportion(dat[,2:p])
    nameSpecies=colnames(propdat)
    #Store(dat)
    
    nbAllSpecies=ncol(propdat)
    
    t1=Sys.time()           

    if (RunHAC == TRUE) {

      # HAC
      print("######## SPECIES EXPLORATION METHOD 1: HAC ########")
      # Transposing data
      table_var=table_variables(propdat)
  
      # HAC
      print("cluster...")
      cah_var=hcluster(table_var, method="euclidean", link="ward")
  
  
      # Select the number of clusters by scree-test
      inerties.vector=cah_var$height[order(cah_var$height,decreasing=T)]
      nb.finalclusters=which(scree(inerties.vector)[,"epsilon"]<0)[1]
  
  
      # Dendogram cutting at the selected level
      cah_cluster_var=cutree(cah_var,k=nb.finalclusters)
  
      png(paste(analysisName,"HAC_Dendogram_Step1.png",sep="_"), width = 1200, height = 800)
      plclust(cah_var,labels=F,hang=-1,ann=F)
      title(main="HAC dendogram",xlab="Species",ylab="Height")
      rect.hclust(cah_var, k=nb.finalclusters)
      dev.off()
        
      temp=select_species(dat[,2:p],cah_cluster_var)
      namesResidualSpecies=nameSpecies[which(cah_cluster_var==temp[[2]])] #list of residual species 
  
      fait=FALSE
      nb_cut=1
      while (( fait == FALSE) && (nb_cut < (p-nb.finalclusters))) {
        # cutting below
        print(paste("-----------   nb_cut =",nb_cut))
        cah_cluster_var_step=cutree(cah_var,k=(nb.finalclusters+nb_cut))
        # testing residual species
        print(paste("------------- Residual species cluster(s) ",unique(cah_cluster_var_step[namesResidualSpecies])))
        if (length(unique(cah_cluster_var_step[namesResidualSpecies]))==1) {
          print(paste("-------------  No residual cut -----"))
          nb_cut = nb_cut+1 #cutting below
        }else{
          print("-------------  Residual cut -----")
          nbSpeciesClusters=table(cah_cluster_var_step[namesResidualSpecies])
          #testing if a species is alone in a group
          if (sort(nbSpeciesClusters)[1]>1) { #if not alone
            print("------- I stop and have a beer ------")
            fait = TRUE # then I stop
          }else{
            print("------ Updating residual species -----")
            nb_cut = nb_cut+1;  # if alone then cutting below and updating nameResidualSpecies to start again
            numGroupSpeciesAlone = as.numeric(names(sort(nbSpeciesClusters)[1]))
            namesSpeciesAlone = names(cah_cluster_var_step)[which(cah_cluster_var_step==numGroupSpeciesAlone)]
            namesResidualSpecies = namesResidualSpecies[ - which(namesResidualSpecies==namesSpeciesAlone)]
            print(paste("---- Adding new species ---",namesSpeciesAlone))
          }
        }
      } # end of while
  
  
      # Dendogram of the first cut in the residual species cluster
      png(paste(analysisName,"HAC_Dendogram_Step1_ResidualSpecies.png",sep="_"), width = 1200, height = 800)
      plclust(cah_var,labels=F,hang=-1,ann=F)
      title(main="HAC dendogram - Step",xlab="Species",ylab="Height")
      rect.hclust(cah_var, k=(nb.finalclusters+nb_cut))
      dev.off()
         
         
      # Selection of main species
      nomespsel=setdiff(nameSpecies,namesResidualSpecies)
      cat("main species : ",nomespsel,"\n")
  
  
      # Return the dataset retaining only the main species
      nbMainSpeciesHAC=length(nomespsel)
      NamesMainSpeciesHAC=nomespsel
      PropMainSpeciesHAC=nbMainSpeciesHAC/nbAllSpecies*100   
      
      if(DiagFlag==TRUE) {
        datSpeciesWithoutProp=building_tab_pca(dat[,2:p],nomespsel)
        pourcentCatchMainSpeciesHAC=apply(datSpeciesWithoutProp,1,sum)/apply(dat[,2:p],1,sum)*100
        medianPourcentCatchMainSpeciesHAC=median(pourcentCatchMainSpeciesHAC)
      }

  print(Sys.time()-t1)

  }else{ NamesMainSpeciesHAC=NA; nbMainSpeciesHAC=as.numeric(NA); medianPourcentCatchMainSpeciesHAC=as.numeric(NA); PropMainSpeciesHAC=NA }
 
    # TOTALE
    
    print("######## SPECIES EXPLORATION METHOD 2: 'TOTAL' ########")
    
    # Total quantity caught by species
    sumcol=numeric(length=p-1)
    for(i in 2:p){
      sumcol[i-1]=sum(dat[,i], na.rm=TRUE)
    }
    names(sumcol)=names(dat)[-1]

    # Percent of each species in the total catch
    propesp=sumcol/sum(sumcol,na.rm=T)*100
    # Columns number of each species by decreasing order of capture
    numesp=order(propesp,decreasing=T)
    # Percent of each species in the total catch by cumulated decreasing order
    propesp=cumsum(propesp[order(propesp,decreasing=T)])

    # We are taking all species until having at least seuil% of total catch

    nbMainSpeciesTotal=numeric()
    medianPourcentCatchMainSpeciesTotal=numeric()

    for(seuil in seq(5,100,5)){
      cat("seuil:",seuil,"\n")
      pourcent=which(propesp<=seuil)
      # We are taking the name of selected species
      espsel=numesp[1:(length(pourcent)+1)]
      nomespsel=nameSpecies[espsel]
#      nomespsel=names(pourcent)
      nbMainSpeciesTotal[seuil/5]=length(nomespsel)

      if(DiagFlag==TRUE) {
        # We are bulding the table with main species and aggregated other species
        datSpeciesWithoutProp=building_tab_pca(dat[,2:p],nomespsel)
        if(length(nomespsel)==1){
          vectorNul=rep(0,n)
          datSpeciesWithoutProp=cbind(datSpeciesWithoutProp,vectorNul)
        }
        pourcentCatchMainSpeciesTotal=apply(datSpeciesWithoutProp,1,sum, na.rm=TRUE)/apply(dat[,2:p],1,sum, na.rm=TRUE)*100
        medianPourcentCatchMainSpeciesTotal[seuil/5]=median(pourcentCatchMainSpeciesTotal)
      }
    }
    nbMainSpeciesTotal=c(0,nbMainSpeciesTotal)
    nbMainSpeciesTotal[length(nbMainSpeciesTotal)]=p-1
    NamesMainSpeciesTotal=nomespsel[1:nbMainSpeciesTotal[length(nbMainSpeciesTotal)-1]]
    PropMainSpeciesTotal=nbMainSpeciesTotal[length(nbMainSpeciesTotal)-1]/nbAllSpecies*100
    
    if (DiagFlag) medianPourcentCatchMainSpeciesTotal=c(0,medianPourcentCatchMainSpeciesTotal)

    print(Sys.time()-t1)



    # LOGEVENT
    
    print("######## SPECIES EXPLORATION METHOD 3: 'LOGEVENT' ########")
    
    nbMainSpeciesLogevent=numeric()
    medianPourcentCatchMainSpeciesLogevent=numeric()

    for(seuil in seq(5,100,5)){
      cat("seuil:",seuil,"\n")
      nomespsel=character()
      for (i in nameSpecies) if (!is.na(any(propdat[,i]>=seuil)) && any(propdat[,i]>=seuil)) nomespsel <- c(nomespsel,i)
      nbMainSpeciesLogevent[seuil/5]=length(nomespsel)
   
      # We are bulding the table with main species and aggregated other species
      if(DiagFlag==TRUE) {
        datSpeciesWithoutProp=building_tab_pca(dat[,2:p],nomespsel)
        if(length(nomespsel)==1){
          vectorNul=rep(0,n)
          datSpeciesWithoutProp=cbind(datSpeciesWithoutProp,vectorNul)
        }
        pourcentCatchMainSpeciesLogevent=apply(datSpeciesWithoutProp,1,sum)/apply(dat[,2:p],1,sum)*100
        medianPourcentCatchMainSpeciesLogevent[seuil/5]=median(pourcentCatchMainSpeciesLogevent)
      }
    }
    nbMainSpeciesLogevent=c(p-1,nbMainSpeciesLogevent)
    NamesMainSpeciesLogevent=nomespsel
    PropMainSpeciesLogevent=nbMainSpeciesLogevent[length(nbMainSpeciesLogevent)]/nbAllSpecies*100
    
    
    if(DiagFlag) medianPourcentCatchMainSpeciesLogevent=c(100,medianPourcentCatchMainSpeciesLogevent)

    print(Sys.time()-t1)

    # GRAPHICS
    # Number of main species
    png(paste(analysisName,"Number of main species.png",sep="_"), width = 1200, height = 800)
    plot(seq(0,100,5),nbMainSpeciesTotal,type='l',col="blue",lwd=2, main="Number of main species selected depending of the threshold", xlab="Threshold (%)",ylab="Number of species")
    lines(seq(0,100,5),nbMainSpeciesLogevent,col="green",lwd=2)
    if(!is.na(nbMainSpeciesHAC)) abline(nbMainSpeciesHAC,0, col="red",lwd=2)
    mtext(paste(p-1," Species"),col='darkblue')
    legend(70, p*0.9, c( "HAC", "Total", "Logevent"),lwd=2,col=c("red", "blue", "green"))
    dev.off()

    # Median percentage of catch represented by main species by logevent
    if(DiagFlag){
      png(paste(analysisName,"Median percentage of catch represented by main species by logevent.png",sep="_"), width = 1200, height = 800)
      plot(seq(0,100,5),medianPourcentCatchMainSpeciesTotal,type='l',col="blue",lwd=2, main="Median percentage of catch represented by main species by logevent depending of the threshold", xlab="Threshold (%)",ylab="Median percentage of catch represented by main species by logevent")
      lines(seq(0,100,5),medianPourcentCatchMainSpeciesLogevent,col="green",lwd=2)
      if (RunHAC==TRUE) abline(medianPourcentCatchMainSpeciesHAC,0, col="red",lwd=2)
      mtext(paste(p-1," Species"),col='darkblue')
      if (RunHAC==TRUE) legend(70, 40, c("HAC", "Total", "Logevent"),lwd=2,col=c("red", "blue", "green"))
      if (RunHAC==FALSE) legend(70, 40, c("Total", "Logevent"),lwd=2,col=c("blue", "green"))
      dev.off()
    }
    
    ListSpecies=sort(unique(c(NamesMainSpeciesHAC,NamesMainSpeciesTotal,NamesMainSpeciesLogevent)))

    # Proportion of the total catch represented by the species in ListSpecies
    catchListSpecies=sumcol[ListSpecies]
    propCatchListSpecies=sum(catchListSpecies)/sum(sumcol)*100
    

    if(DiagFlag==FALSE) { 
      explo_species = list(nbAllSpecies=nbAllSpecies,
                            PropMainSpeciesHAC=PropMainSpeciesHAC,
                            PropMainSpeciesTotal=PropMainSpeciesTotal,
                            PropMainSpeciesLogevent=PropMainSpeciesLogevent,
                            nbMainSpeciesHAC=nbMainSpeciesHAC, 
                            nbMainSpeciesTotal=nbMainSpeciesTotal, 
                            nbMainSpeciesLogevent=nbMainSpeciesLogevent,
                            NamesMainSpeciesHAC=sort(NamesMainSpeciesHAC), 
                            NamesMainSpeciesTotalAlphabetical=sort(NamesMainSpeciesTotal),                                             
                            NamesMainSpeciesTotalByImportance=NamesMainSpeciesTotal,
                            NamesMainSpeciesLogevent=sort(NamesMainSpeciesLogevent),
                            NamesMainSpeciesAll=ListSpecies,
                            propCatchMainSpeciesAll=propCatchListSpecies) 
    }else{         
      explo_species = list(nbAllSpecies=nbAllSpecies,
                            PropMainSpeciesHAC=PropMainSpeciesHAC,
                            PropMainSpeciesTotal=PropMainSpeciesTotal,
                            PropMainSpeciesLogevent=PropMainSpeciesLogevent,
                            nbMainSpeciesHAC=nbMainSpeciesHAC, 
                            nbMainSpeciesTotal=nbMainSpeciesTotal, 
                            nbMainSpeciesLogevent=nbMainSpeciesLogevent,
                            NamesMainSpeciesHAC=sort(NamesMainSpeciesHAC), 
                            NamesMainSpeciesTotalAlphabetical=sort(NamesMainSpeciesTotal),                                             
                            NamesMainSpeciesTotalByImportance=NamesMainSpeciesTotal,
                            NamesMainSpeciesLogevent=sort(NamesMainSpeciesLogevent),
                            NamesMainSpeciesAll=ListSpecies,
                            medianPourcentCatchMainSpeciesHAC=median(pourcentCatchMainSpeciesHAC),
                            medianPourcentCatchMainSpeciesTotal=medianPourcentCatchMainSpeciesTotal,
                            medianPourcentCatchMainSpeciesLogevent=medianPourcentCatchMainSpeciesLogevent,
                            propCatchMainSpeciesAll=propCatchListSpecies)
    }    

    return(explo_species)
    
}