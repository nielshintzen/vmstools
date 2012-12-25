##########################################################################
#    multivariateClassificationUtils.R  CONTAINS ALL USEFUL FUNCTIONS    #
#    FOR THE MULTIVARIATE CLASSIFICATION                                 #
##########################################################################


require(FactoMineR)   # function PCA
require(cluster)      # functions pam & clara
require(SOAR)         # function Store
require(amap)         # function hcluster
require(MASS)         # function lda
require(mda)          # function fda



# Transform quantities to percentage values (between 0 and 100) 
# of each species in the logevent total catch
transformation_proportion=function(tab){
  res=as.matrix(tab)
  n=nrow(tab)
  p=ncol(tab)
  for (i in 1:n){
    sommeligne=sum(res[i,], na.rm=TRUE)
    if(sommeligne==0){
      res[i,]=rep(0,p)
    }else{
      res[i,]=res[i,]*(100/sommeligne)
    }
  }
  return(res)
}


# Transpose the dataset (change variables into individuals)
table_variables=function(data){
  n=nrow(data)
  res1=t(as.matrix(data[1:round(n/2),]))
  res2=t(as.matrix(data[(round(n/2)+1):n,]))
  res=cbind(res1,res2)
  row.names(res)=colnames(data)
  colnames(res)=row.names(data)
  return(res)
}


# Implementation of "scree-test"
scree=function(eig){
  n=length(eig)
  delta=numeric(n)
  epsilon=numeric(n)
  delta[2]=eig[2]-eig[1]
  for (i in 3:n){
    delta[i]=eig[i]-eig[i-1]
    epsilon[i]=delta[i]-delta[i-1]
  }
  data=matrix(0,nrow=n,ncol=3)
  data=cbind(valeurs_propres=eig, delta=delta, epsilon=epsilon)
  return(data)
}


# Remove the cluster with the smallest mean of capture
select_species=function(data,groupes_cah){
  nb.classes=length(levels(as.factor(groupes_cah)))
  moyennes=numeric(nb.classes)
  for(i in 1:nb.classes){
    namegp=names(which(groupes_cah==i))
    effgp=length(which(groupes_cah==i))
    moyennes[i]=sum(data[namegp],na.rm=TRUE)/effgp
  }
  indice.autre=which(moyennes == min(moyennes,na.rm=TRUE)) 
  noms=names(which(groupes_cah!=indice.autre))
  return(list(noms,indice.autre))
}


# Build the table with the main species
building_tab_pca=function(data,especes){
  p=ncol(data)
  noms=colnames(data)
  ind_princ=which(is.element(noms,especes))
  princ=data[,ind_princ]
  return(princ)
}


# Compute the test-value for each species by cluster
test.values=function(groupes,data){
                                                                         
  n=nrow(data)
  p=ncol(data)
  noms_var=colnames(data)
  nb_groupes=length(levels(as.factor(groupes)))
  noms_groupes=character(nb_groupes)

  stats_globales=matrix(0,nr=p,ncol=2)
  row.names(stats_globales)=noms_var
  colnames(stats_globales)=c("mean","variance")
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
      mu=stats_globales[noms_var[i],"mean"]
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
    noms_groupes[j]=paste("Cluster",j,sep=" ")
  }
  colnames(res)=noms_groupes
  return(res)
}


# Determine the species with a test-value > 1.96 by cluster
targetspecies=function(resval){
  p=nrow(resval)
  nbgp=ncol(resval)
  
  tabnumespcib=data.frame()
  tabnomespcib=data.frame()
   
  for(i in 1:nbgp){
    # qnorm(0.975,mean=0,sd=1)=1.96     (P(resval>1.96)=0.025)
    numespcib=which(resval[,i]>1.96)   
    numespcibdec=numespcib[order(resval[numespcib,i],decreasing=TRUE)]           
    nomespcib=names(numespcibdec)
       
    nbespgpcib=length(numespcib)

    if(nbespgpcib>0){
      for (j in 1:nbespgpcib){
        tabnumespcib[i,j]=numespcibdec[j]
        tabnomespcib[i,j]=nomespcib[j]
      }                                                       
    }else{
        tabnumespcib[i,]=NA
        tabnomespcib[i,]=NA
    }   
  }
  tabnumespcib=as.matrix(tabnumespcib)
  tabnomespcib=as.matrix(tabnomespcib)
  return(list(tabnumespcib=tabnumespcib,tabnomespcib=tabnomespcib))
}


