#' Finding metiers for a new reduced EFLALO dataset from a reduced EFLALO
#' dataset for which the metiers are known.
#' 
#' This function permits to determine the metiers of new logevents belonging to
#' a new reduced EFLALO dataset, from a reduced EFLALO dataset for which the
#' metiers are known. For example, this function permits to determine the
#' metiers of logevents realized in 2008, basing on metiers of logevents
#' realized in 2007. It works with a discriminante analysis ran with the
#' function fda. In the first time, the function fda calibrates a model basing
#' on the dataset for which the metiers are known. In the second time, the
#' function predict is used to predict the metiers of new logevents, thanks to
#' the previous model.
#' 
#' 
#' @param learningData A numerical matrix with logevents as lines and species
#' as columns, with percentage values (between 0 and 100) of each species in
#' the logevent catches. Logevent ID (LE_ID) should be as row names. Typically,
#' this table will be produced from a eflalo dataset at the step 1 of the
#' metier analysis, using the function extractTableMainSpecies()
#' @param clustersAffectation An integer vector of length n, the number of
#' logevents in learningData, giving for each logevent the number of the
#' cluster to which it belongs. This vector will be produced at the step 3 of
#' the metier analysis, using the function getMetierClusters().
#' @param newData A numerical matrix with logevents as lines and species as
#' columns, with raw values of each species in the logevent catches.  Logevent
#' ID (LE_ID) should be as row names.
#' @return A data frame giving the logevent ID (LE_ID) (first column) and the
#' number of the cluster to which it belongs (second column), for all logevents
#' in newData. The cluster number "0" corresponds to the unknown metier.
#' @note A number of libraries are initially called for the whole metier
#' analyses and must be installed :
#' (FactoMineR),(cluster),(SOAR),(amap),(MASS),(mda)
#' @author Nicolas Deporte, Sebastien Demaneche, Stephanie Mahevas (IFREMER,
#' France), Clara Ulrich, Francois Bastardie (DTU Aqua, Denmark)
#' @seealso \code{\link{selectMainSpecies}},
#' \code{\link{extractTableMainSpecies}}, \code{\link{getTableAfterPCA}},
#' \code{\link{getMetierClusters}}
#' @references Development of tools for logbook and VMS data analysis. Studies
#' for carrying out the common fisheries policy No MARE/2008/10 Lot 2
#' @examples
#' 
#' 
#'   \dontrun{
#' 
#'   data(eflalo)
#' 
#'   # Format
#'   newEflalo <- formatEflalo(eflalo)
#'   newEflalo <- newEflalo[newEflalo$LE_GEAR=="OTB",]
#'   
#'   # Note that output plots will be sent to getwd()
#'   analysisName <- "metier_analysis_OTB"
#' 
#'   # Prepare datPred
#'   datPred <- newEflalo[,c("LE_ID",grep("EURO",names(newEflalo),value=TRUE))]
#'   datPred[is.na(datPred)] <- 0
#'   names(datPred)[-1] <- unlist(lapply(strsplit(names(datPred[,-1]),"_"),
#'                                function(x) x[[3]]))
#'   le_id_datPred <- datPred[,"LE_ID"]
#'   datPred <- datPred[,-1]
#'   datPred <- as.matrix(datPred)
#'   rownames(datPred) <- le_id_datPred
#' 
#'   if (methMetier=="hac")     clustersAffectation <- Step3$clusters
#'   if (methMetier=="clara")   clustersAffectation <- Step3$clusters$clustering
#'   if (methMetier=="kmeans")  clustersAffectation <- Step3$clusters$cluster
#' 
#'   # Predict the metiers of the new logbook events of newData from the metiers
#'   # found in Step 3
#'   metierPred <- predictMetier(learningData=Step1,
#'                   clustersAffectation=clustersAffectation,newData=datPred)
#'     
#'   }
#' 
#' 
#' @export predictMetier
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
