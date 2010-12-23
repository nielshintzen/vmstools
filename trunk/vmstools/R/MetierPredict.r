####################################################
#   DISCRIMINANTE ANALYSIS FOR METIER PREDICTION   #
####################################################

metierPredict=function(learningData,clustersAffectation,newData){
#,clustersCheck){

  print("#--------- Predict Metier --------#")
  #library(MASS)
  library(klaR)   # rda
  #library(mda)   # fda

  le_id=rownames(learningData)
  nbSpeciesLearningData=ncol(learningData)   # Number of species of learningData
  namelearningDataSpecies=colnames(learningData)
  #learningData=as.matrix(learningData,ncol=nbSpeciesLearningData,nrow=length(le_id))

  # Discriminante analysis on the learning dataset
  #learning.lda=lda(learningData,as.factor(clustersAffectation))
  
  learningData=as.data.frame(cbind(learningData,as.factor(clustersAffectation)))
  colnames(learningData)=c(colnames(learningData[1:nbSpeciesLearningData]),"clust")
  print("------------ Learning ------------")
  learning=rda(clust~.,data=learningData,gamma=1,lambda=1)
  #learning=fda(clust~.,data=learningData,method=mars)
  
  #Moving  logevents of newData characterized by catch of species not caught in datSpecies
  #in an undefined Metier

  nbSpecNewData=ncol(newData)   # Number of species of newData
  nbLogNewData=nrow(newData)   #Number of logevents of newData
  nameNewDataSpecies=colnames(newData)
  #list of species in NewData not present in learningDataSpecies
  newSpecies=nameNewDataSpecies[is.na(match(nameNewDataSpecies,namelearningDataSpecies))]
  
  #logevents with positive catch for at least one of newSpecies
  nameLogEventNewData=rownames(newData)
  temp=(newData[,newSpecies]>0)
  LE_ID_NewData_unknownMetier=nameLogEventNewData[apply(temp,1,sum)>0]
  notPredictedClusters=data.frame(LE_ID_NewData_unknownMetier,
        as.factor(rep(0,length(LE_ID_NewData_unknownMetier))))
  dimnames(notPredictedClusters)[[2]]=c("LE_ID","Class")
  #format newData to get all species in datSpecies an removing species of newData not present in datSpecies
  newDataFormat0=newData[-(which(is.element(nameLogEventNewData,LE_ID_NewData_unknownMetier))),
  -which(is.element(nameNewDataSpecies,newSpecies))]
  
  newDataFormat=merge(newDataFormat0,learningData,all=TRUE)
  newDataFormat=newDataFormat[which(is.element(rownames(newDataFormat0),nameLogEventNewData)),]
  newDataFormat[is.na(newDataFormat)]=0


  #newDataFormat=newDataFormat[,2:p] #remove LE_ID
  newOrder=match(namelearningDataSpecies,colnames(newDataFormat))
#  newDataFormatMatrix = as.matrix(newDataFormat[,newOrder],
#            ncol=nbSpeciesLearningData , nrow=nbLogNewData) #Order col as in learningData
  newDataFormatMatrix=as.data.frame(newDataFormat[,newOrder],
            ncol=nbSpeciesLearningData , nrow=nbLogNewData) #Order col as in learningData
  
  #predict clusters
  #result.lda=predict(learning.lda,newDataFormatMatrix)
  print("------------ Predicting ------------")
  result=predict(learning,newDataFormatMatrix)
  #predictedClusters=data.frame(rownames(newDataFormat0),result.lda$class)
  predictedClusters=data.frame(rownames(newDataFormat0),result$class)
  dimnames(predictedClusters)[[2]]=c("LE_ID","Class")

  clustersForAllLogEvents=rbind(predictedClusters,notPredictedClusters)

  return(clustersForAllLogEvents=clustersForAllLogEvents)
}