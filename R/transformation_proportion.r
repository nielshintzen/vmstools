#################################################################
# Transform quantities to percentage values (between 0 and 100) #
#        of each species in the logevent total catch            #
#################################################################

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