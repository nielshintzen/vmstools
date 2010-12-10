####################################################
#   DISCRIMINANTE ANALYSIS FOR METIER PREDICTION   #
####################################################

metierPredict=function(learningData,clustersAffectation,newData){
#,clustersCheck){

  library(MASS)

  p=ncol(learningData)   # Number of species of learningData +1

  # Discriminante analysis on the learning dataset
  learning.lda=lda(learningData[,2:p],as.factor(clustersAffectation))

  #Moving  logevents of newData characterized by catch of species not caught in datSpecies
  #in an undefined Metier

   namelearningDataSpecies=colnames(learningData)[2:p]
   q=ncol(newData)   # Number of species of newData  +1
   l=nrow(newData)   #Number of logevents of newData
   nameNewDataSpecies = colnames(newData)[2:q]
   #list of species in NewData not present in learningDataSpecies
   newSpecies=nameNewDataSpecies[is.na(match(nameNewDataSpecies,namelearningDataSpecies))]
   
   #logevents with positive catch for at least one of newSpecies
   temp=(newData[,newSpecies]>0)
   LE_ID_NewData_unknownMetier=newData$LE_ID[apply(temp,1,sum)>0]
   #number of unknown cluster = 0 
   notPredictedClusters = cbind(LE_ID_NewData_unknownMetier,rep(0,length(LE_ID_NewData_unknownMetier)))
   
   #format newData to get all species in datSpecies an removing species of newData not present in datSpecies
   newDataFormat0=newData[-(which(is.element(newData$LE_ID,LE_ID_NewData_unknownMetier))),
    -which(is.element(colnames(newData),newSpecies))]

   #NbLogToBePredicted=newDataFormat0[[1]]
   #RemainingNewDataSpecies=colnames(newDataFormat0)[-1]


   newDataFormat=merge(newDataFormat0,learningData,all=TRUE)
   newDataFormat=newDataFormat[which(is.element(newDataFormat0$LE_ID,newDataFormat$LE_ID)),]
   newDataFormat[is.na(newDataFormat)]=0


  newDataFormat=newDataFormat[,2:p] #remove LE_ID
  newOrder=match(colnames(learningData)[-1],colnames(newDataFormat))
  newDataFormatMatrix = as.matrix(newDataFormat[,newOrder]) #Order col as in learningData
  
  #predict clusters
  result.lda=predict(learning.lda,newDataFormatMatrix)
  predictedClusters=cbind(newDataFormat0$LE_ID,result.lda$class)


  clustersForAllLogEvents = rbind(predictedClusters,notPredictedClusters)
  
  return(clustersForAllLogEvents=clustersForAllLogEvents)
}