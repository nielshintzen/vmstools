#############################################################
# Determine the species with a test-value > 1.96 by cluster #
#############################################################

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

