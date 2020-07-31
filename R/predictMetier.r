predictMetier=function(learningData=Step1,clustersAffectation=clust2007,newData=datPred){

    # Select only the columns of newData corresponding to species of learningData
    namesSpecies=colnames(learningData)
    le_id_newData=rownames(newData)
    newData=newData[,namesSpecies]
    # Transform quantities of newData to percentages of each logevent totale catch
    newData=transformation_proportion(newData)

    # Select the logevents without catch for the selected species
    nullCatch=which(apply(newData,1,sum)==0)
    newDataWithoutCatch=newData[nullCatch,]
    if(length(nullCatch)==1){
      le_id_newDataWithoutCatch=names(nullCatch)
    }else{
      le_id_newDataWithoutCatch=rownames(newDataWithoutCatch)
    }

    # Select the logevents with catch for the selected species
    positiveCatch=setdiff(1:nrow(newData),nullCatch)
    newDataWithCatch=newData[positiveCatch,]
    
    # Prepare learningData for the discriminant analysis
    le_id_learningData=rownames(learningData)
    nbSpeciesLearningData=ncol(learningData)   # Number of species of learningData
    nameLearningDataSpecies=colnames(learningData)
    learningData=as.data.frame(cbind(learningData,as.factor(clustersAffectation)))
    colnames(learningData)=c(colnames(learningData[1:nbSpeciesLearningData]),"clust")

    # Calibrate the model with learningData
    print("------------ Learning ------------")
    learning=fda(clust~.,data=learningData)
    
    # Predict the metier of each logevent in newDataWithCatch
    print("------------ Predicting ------------")
    result=predict(learning,newdata=newDataWithCatch)
    predictedClusters=data.frame(rownames(newDataWithCatch),result)
    dimnames(predictedClusters)[[2]]=c("LE_ID","Class")
    
    # Give the metier "0" (unknown metier) for the logevents in newDataWithoutCatch
    if(length(le_id_newDataWithoutCatch)!=0){
      notPredictedClusters=data.frame(le_id_newDataWithoutCatch,
          as.factor(rep(0,length(le_id_newDataWithoutCatch))))
      dimnames(notPredictedClusters)[[2]]=c("LE_ID","Class")
    }
    
    # Create a table linking "LE_ID" and metier for each logevent in newData
    if(length(le_id_newDataWithoutCatch)!=0){
      clustersForAllLogevents=rbind(predictedClusters,notPredictedClusters)
    }else{
      clustersForAllLogevents=predictedClusters
    }
    
    return(clustersForAllLogevents=clustersForAllLogevents)
}