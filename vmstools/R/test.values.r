######################################################
# Compute the test-value for each species by cluster #
######################################################

test.values=function(groupes,data){

  n=nrow(data)
  p=ncol(data)
  noms_var=colnames(data)
  nb_groupes=length(levels(as.factor(groupes)))
  noms_groupes=character(nb_groupes)

  stats_globales=matrix(0,nrow=p,ncol=2)
  row.names(stats_globales)=noms_var
  colnames(stats_globales)=c("mean","variance")
  for (i in 1:p){
    stats_globales[i,1]=mean(data[,noms_var[i]])
    stats_globales[i,2]=var(data[,noms_var[i]])
  }

  res=matrix(0,nrow=p,ncol=nb_groupes)
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
