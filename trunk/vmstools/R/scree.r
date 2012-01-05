###################################
# Implementation of "scree-test"  #
###################################

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
