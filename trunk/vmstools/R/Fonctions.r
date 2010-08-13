###############################################################################
#    Fonctions.r  contains all useful functions for the function "classif"    #
###############################################################################


library(FactoMineR)
library(ade4)
library(cluster)
require(graphics)
require(stats)
library(SOAR)
library(amap)
#library(spdep)   (function ssw)
library(MASS)


# Transform quantities to percents fished in the logevent

transformation_proportion=function(tab){
  res=as.matrix(tab)
  n=length(tab[,1])
  p=length(tab[1,])
  for (i in 1:n){
    sommeligne=sum(res[i,])
    if(sommeligne==0){
      res[i,]=0
    }else{
      res[i,]=res[i,]*(100/sommeligne)
    }
  }
  return(data.frame(res))
}

#transformation_proportion=function(tab){
#  res=as.matrix(tab)
#  res=apply(tab,1, function(x) { if(sum(x)==0){ 0 }else{ x/sum(x)*100 } })
#  res=t(res)
#  return(data.frame(res))
#}


# Transposing data (changing variables into individuals)

table_variables=function(data){
  n=length(data[,1])
  p=length(data[1,])
  res1=t(as.matrix(data[1:round(n/2),]))
  res2=t(as.matrix(data[(round(n/2)+1):n,]))
  res=cbind(res1,res2)
  res=data.frame(res)
  row.names(res)=colnames(data)
  colnames(res)=row.names(data)
  return(res)
}


# Scree-test

scree=function(eig){
  n=length(eig)
  delta=numeric(n)
  epsilon=numeric(n)
  delta[2]=eig[2]-eig[1]
  for (i in 3:n){
    delta[i]=eig[i]-eig[i-1]
    epsilon[i]=delta[i]-delta[i-1]
  }

  data=data.frame(valeurs_propres=eig, delta=delta, epsilon=epsilon)
  return(data)
}


# Removing the cluster with the smaller mean of capture

select_species=function(data,groupes_cah){
  nb.classes=length(levels(as.factor(groupes_cah)))
  moyennes=numeric(nb.classes)
  for(i in 1:nb.classes){
    namegp=names(which(groupes_cah==i))
    effgp=length(which(groupes_cah==i))
    moyennes[i]=sum(data[namegp])/effgp
  }
  indice.autre=which.min(moyennes)
  noms=names(which(groupes_cah!=indice.autre))
  return(noms)
}


# Building the table with main species and agregated OtherSpecies

building_tab_pca=function(data,especes){
  p=ncol(data)
  noms=names(data)
  ind_princ=which(is.element(noms,especes))
  ind_autres=setdiff(1:p,ind_princ)
  princ=data[,ind_princ]
  autres=apply(data[,ind_autres],1,sum)
  #return(data.frame(princ,OTHERS=autres))
  return(data.frame(princ))
}


# Function computing the test-values

test.values=function(groupes,data){
                                                                         
  n=nrow(data)
  p=ncol(data)
  noms_var=colnames(data)
  nb_groupes=length(levels(as.factor(groupes)))

  noms_groupes=character(nb_groupes)

  stats_globales=matrix(0,nr=p,ncol=2)
  row.names(stats_globales)=noms_var
  colnames(stats_globales)=c("moyenne","variance")
  for (i in 1:p){
    stats_globales[i,1]=mean(data[,noms_var[i]])
    stats_globales[i,2]=var(data[,noms_var[i]])
  }

  res=matrix(0,nr=p,ncol=nb_groupes)
  row.names(res)=noms_var

  for (j in 1:nb_groupes){
    groupe=which(groupes==j)                        
    n_k=length(groupe)

    for (i in 1:p){
      mu_k=mean(data[groupe,noms_var[i]])
      mu=stats_globales[noms_var[i],"moyenne"]
      V=stats_globales[noms_var[i],"variance"]
      V_mu_k=(n-n_k)*V/(n_k*(n-1))
      
      if(V_mu_k==0){
        Valeur_test=0
      }else{
        Valeur_test=(mu_k-mu)/sqrt(V_mu_k)
      }
                                                
      res[i,j]=Valeur_test
      rm(Valeur_test)
    }
    rm(groupe)
    noms_groupes[j]=paste("groupe",j,sep="_")
  }
  colnames(res)=noms_groupes
  return(data.frame(res))
}


# Fonction determining the target species

targetspecies=function(resval){
  p=length(resval[,1])
  nbgp=length(resval[1,])
  
  tabnumespcib=data.frame()
  tabnomespcib=data.frame()
  tabnumespsignif=data.frame()
  tabnomespsignif=data.frame()
  tabnumesppascib=data.frame()
  tabnomesppascib=data.frame()
    
  for(i in 1:nbgp){
    # qnorm(0.975,mean=0,sd=1)=1.96
    numespcib=which(resval[,i]>1.96)   
    numespcibdec=numespcib[order(resval[numespcib,i],decreasing=T)]           
    nomespcib=row.names(resval[numespcibdec,])
    
    numespsignif=which(abs(resval[,i])>1.96)
    nomespsignif=row.names(resval[numespsignif,])
    
    numesppascib=setdiff(numespsignif,numespcib)
    numesppascibdec=numesppascib[order(resval[numesppascib,i],decreasing=F)]
    nomesppascib=row.names(resval[numesppascibdec,])
    
    
    nbespgpcib=length(numespcib)
    nbespgpsignif=length(numespsignif)
    nbespgppascib=length(nomesppascib)

    if(nbespgpcib>0){
      for (j in 1:nbespgpcib){
        tabnumespcib[i,j]=numespcibdec[j]
        tabnomespcib[i,j]=nomespcib[j]
      }                                                       
    }else{
        tabnumespcib[i,]=NA
        tabnomespcib[i,]=NA
    }
    
    if(nbespgpsignif>0){
      for (j in 1:nbespgpsignif){
        tabnumespsignif[i,j]=numespsignif[j]
        tabnomespsignif[i,j]=nomespsignif[j]
      }
    }else{
        tabnumespsignif[i,]=NA
        tabnomespsignif[i,]=NA
    }
    
    if(nbespgppascib>0){
      for (j in 1:nbespgppascib){
        tabnumesppascib[i,j]=numesppascibdec[j]
        tabnomesppascib[i,j]=nomesppascib[j]
      }
    }else{
        tabnumesppascib[i,]=NA
        tabnomesppascib[i,]=NA
    }          
    
  }
  return(list(tabnumespcib=tabnumespcib,tabnomespcib=tabnomespcib))
}


# Calcule of within variance of a cluster (square distance between each row of the cluster and the cluster's center of gravity)
# we calculate the distance row by row
withinVar=function(oneRowOfCluster,centerOfGravityClusti){
  CoGWithOneRow=numeric()
  CoGWithOneRow=rbind(centerOfGravityClusti, oneRowOfCluster)
  distRowCoG=dist(CoGWithOneRow)^2
  return(distRowCoG)
}




# Calculation of the selection's threshold of main species for method "logevent"
#
#seuillogevent=function(propdat){
#  p=length(propdat[1,])   # number of species
#  hcounts=data.frame()
#  fc=data.frame()
#  diffc=data.frame()
#  deltamax=numeric()
#  indmax=numeric()
#  
#  for (i in 1:p){
#    his=hist(propdat[,i], plot=F)$counts
#    if(length(his)==1){
#      his=c(his,rep(0,9))
#    }
#    hcounts=rbind(hcounts,his)   # abscisse : % de sole dans la capture // ordonnée : nombre de logevents pour lesquels la sole représente entre tant et tant de % de la capture totale
#    freqc=cumsum(his/sum(his)*100)  # abscisse :  % de sole dans la capture // ordonnée : % des logevents pour lesquels la sole représente moins de tant de % de la capture totale
#    fc=rbind(fc,freqc)
#    #x=seq(10,100,10)    # abscisses
#    #plot(x,fc1)
#    freqc=c(freqc[1],freqc)
#    diffreqc=diff(freqc)
#    diffc=rbind(diffc,diffreqc)
#    #plot(x,diffc1)
#    deltamax[i]=max(diffreqc)
#    indmax[i]=which(diffreqc==deltamax[i])*10
#  }
#  
#  # pb au niveau du calcul des différences de fréquences cumulées  
#  pb=which(fc[,10]-fc[,9]<0)        
#  ok=setdiff(1:p,pb)  
#  
#  Delta=mean(indmax,na.rm=T)
#  espprinc=which(indmax>Delta)                    
#  Deltamin=min(indmax[espprinc])
#  
#  return(Deltamin)
#}

