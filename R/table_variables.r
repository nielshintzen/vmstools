#############################################################
# Transpose the dataset (change variables into individuals) #
#############################################################

table_variables=function(data){
  n=nrow(data)
  res1=t(as.matrix(data[1:round(n/2),]))
  res2=t(as.matrix(data[(round(n/2)+1):n,]))
  res=cbind(res1,res2)
  row.names(res)=colnames(data)
  colnames(res)=row.names(data)
  return(res)
}
