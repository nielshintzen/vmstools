# from step 7 to step 5

path <- "C:/CLARA/VMSTools_sandbox/"
#path <- "X:/VMSTools_sandbox/"

setwd(path)

#extracting the step 3

#load("DK_OTB2007_3a4/PCA_70/CLARA/Step3.Rdata")
load("All_OTB2007_3a4/PCA_70/CLARA/Step3.Rdata")

# characteristing the main species

#checking the most robust metrics

## I. MINIMUM PERCENTAGE IN THE CLUSTER

#windows(record=TRUE)
#par(mfrow=c(2,2))

#for (metrics in c("1st Qu.","Mean","Median","3rd Qu.")) {
#  mean_ <- Step3$summaryClusters[metrics,,]
#  
#  #plotting the number of main species by percentage threshold
#  threshold <- seq(50,5,by=-5)
#  nb_spp <- matrix(0,ncol=Step3$nbClust,nrow=length(threshold),dimnames=list(threshold,1:Step3$nbClust))
#  for (i in threshold) {
#    spp <- lapply(colnames(mean_),function(x) which(mean_[,x]>=i))
#    nb_spp[as.character(i),] <- sapply(spp, length)
#  }
#  barplot(t(nb_spp),beside=T, ylab="Nb Target Species",xlab="min % threshold in cluster",col=rainbow(Step3$nbClust),main=metrics)
#
#}

# OU : Permet de conserver nb_spp pour chaque metrics.
windows(record=TRUE)
par(mfrow=c(2,2))
threshold <-seq(50,5,by=-5)
nb_spp=array(0,dim=c(length(threshold),Step3$nbClust,4))
dimnames(nb_spp)[[1]]=threshold
dimnames(nb_spp)[[2]]=1:Step3$nbClust
dimnames(nb_spp)[[3]]=c("1st Qu.","Mean","Median","3rd Qu.")

for (metrics in c("1st Qu.","Mean","Median","3rd Qu.")) {
  mean_ <- Step3$summaryClusters[metrics,,]
  #plotting the number of main species by percentage threshold
  for (i in threshold) {
    spp <- lapply(colnames(mean_),function(x) which(mean_[,x]>=i))
    nb_spp[as.character(i),,metrics] <- sapply(spp, length)
  }
  barplot(t(nb_spp[,,metrics]),beside=T, ylab="Nb Target Species",xlab="min % threshold in cluster",col=rainbow(Step3$nbClust),main=metrics)
}



## II. CUMULATED PERCENTAGE IN THE CLUSTER

#for (metrics in c("1st Qu.","Mean","Median","3rd Qu.")) {
#mean_ <- Step3$summaryClusters[metrics,,]
#
##cumulated sum
#threshold <-seq(50,95,by=5) 
#nb_spp <- matrix(0,ncol=Step3$nbClust,nrow=length(threshold),dimnames=list(threshold,1:Step3$nbClust))
#for (i in threshold) {
#spp <- lapply(colnames(mean_),function(x) {
#                                cum_ <- cumsum(sort(mean_[,x],decreasing=TRUE)) 
#                                res <- cum_[seq(length(which(cum_<=i))+1)]
#                                })
#                                
#nb_spp[as.character(i),] <- sapply(spp,length)
#}
#barplot(t(nb_spp),beside=T, ylab="Nb Target Species",xlab="cumulated % in cluster",col=rainbow(Step3$nbClust),main=metrics)
#
#}

# OU : Permet de conserver nb_spp pour chaque metrics.
windows(record=TRUE)
par(mfrow=c(2,2))
threshold <- seq(50,95,by=5)
nb_spp=array(0,dim=c(length(threshold),Step3$nbClust,4))
dimnames(nb_spp)[[1]]=threshold
dimnames(nb_spp)[[2]]=1:Step3$nbClust
dimnames(nb_spp)[[3]]=c("1st Qu.","Mean","Median","3rd Qu.")

for (metrics in c("1st Qu.","Mean","Median","3rd Qu.")) {
  mean_ <- Step3$summaryClusters[metrics,,]
  #cumulated sum 
  for (i in threshold) {
    spp <- lapply(colnames(mean_),function(x) {
                                    cum_ <- cumsum(sort(mean_[,x],decreasing=TRUE)) 
                                    res <- cum_[seq(length(which(cum_<=i))+1)]
                                    })                                
    nb_spp[as.character(i),,metrics] <- sapply(spp,length)
  }
  barplot(t(nb_spp[,,metrics]),beside=T, ylab="Nb Target Species",xlab="cumulated % in cluster",col=rainbow(Step3$nbClust),main=metrics)

}


## III. TEST-VALUE IN THE CLUSTER
windows(record=TRUE)
testValues_=Step3$testValues
threshold <-seq(0,200,by=10)
nb_spp <- matrix(0,ncol=Step3$nbClust,nrow=length(threshold),dimnames=list(threshold,1:Step3$nbClust))
for (i in threshold) {
  spp <- lapply(colnames(testValues_),function(x) which(testValues_[,x]>=i))
  nb_spp[as.character(i),] <- sapply(spp, length)
}
barplot(t(nb_spp),beside=T, ylab="Nb Target Species",xlab="min test-value in cluster",col=rainbow(Step3$nbClust),main="Test-Value")




## IV. % OF LOGEVENT IN THE CLUSTER
windows(record=TRUE)
percLogevents_=Step3$percLogevents
threshold <-seq(0,100,by=5)
nb_spp <- matrix(0,ncol=Step3$nbClust,nrow=length(threshold),dimnames=list(threshold,1:Step3$nbClust))
for (i in threshold) {
  spp <- lapply(rownames(percLogevents_),function(x) which(percLogevents_[x,]>=i))
  nb_spp[as.character(i),] <- sapply(spp, length)
}
barplot(t(nb_spp),beside=T, ylab="Nb Target Species",xlab="min % logevents in cluster",col=rainbow(Step3$nbClust),main="% Logevents")





## V. RESUME
windows(record=TRUE)
par(mfrow=c(3,1))
# par(mfrow=c(2,2))

# % CATCH
#mean_=Step3$summaryClusters["Mean",,]
#thresholdMean=seq(5,50,by=5)
#nb_spp <- matrix(0,ncol=Step3$nbClust,nrow=length(thresholdMean),dimnames=list(thresholdMean,1:Step3$nbClust))
#for (i in thresholdMean) {
#  spp <- lapply(colnames(mean_),function(x) which(mean_[,x]>=i))
#  nb_spp[as.character(i),] <- sapply(spp, length)
#}
#barplot(t(nb_spp),beside=T, ylab="Nb Target Species",xlab="% catch in cluster",col=rainbow(Step3$nbClust),main="% Catch")

# % CUMULATED CATCH
thresholdCumMean=seq(50,95,by=5) 
nb_spp <- matrix(0,ncol=Step3$nbClust,nrow=length(thresholdCumMean),dimnames=list(thresholdCumMean,1:Step3$nbClust))
for (i in thresholdCumMean) {
spp <- lapply(colnames(mean_),function(x) {
                                cum_ <- cumsum(sort(mean_[,x],decreasing=TRUE)) 
                                res <- cum_[seq(length(which(cum_<=i))+1)]
                                })
                                
nb_spp[as.character(i),] <- sapply(spp,length)
}
barplot(t(nb_spp),beside=T, ylab="Number of Target Species",xlab="% Cumulated Catch in cluster",col=rainbow(Step3$nbClust),ylim=c(0,25))

# TEST-VALUES
testValues_=Step3$testValues
thresholdTestValue=seq(0,100,by=10)
nb_spp <- matrix(0,ncol=Step3$nbClust,nrow=length(thresholdTestValue),dimnames=list(thresholdTestValue,1:Step3$nbClust))
for (i in thresholdTestValue) {
  spp <- lapply(colnames(testValues_),function(x) which(testValues_[,x]>=i))
  nb_spp[as.character(i),] <- sapply(spp, length)
}
barplot(t(nb_spp),beside=T, ylab="Number of Target Species",xlab="Test-Value in cluster",col=rainbow(Step3$nbClust),ylim=c(0,25))

# % LOGEVENTS
percLogevents_=Step3$percLogevents
thresholdLogevents <-seq(10,100,by=10)
nb_spp <- matrix(0,ncol=Step3$nbClust,nrow=length(thresholdLogevents),dimnames=list(thresholdLogevents,1:Step3$nbClust))
for (i in thresholdLogevents) {
  spp <- lapply(rownames(percLogevents_),function(x) which(percLogevents_[x,]>=i))
  nb_spp[as.character(i),] <- sapply(spp, length)
}
barplot(t(nb_spp),beside=T, ylab="Number of Target Species",xlab="% Logevents in cluster",col=rainbow(Step3$nbClust),ylim=c(0,25))







# SELECTION OF TARGET SPECIES

# Selection des espèces cibles dans l'ensemble des espèces présentes
#thresholdCatch=75
#thresholdTestValue=1.96
#thresholdLogevents=50
#
#mean_=Step3$summaryClusters["Mean",,]
#testValues_=Step3$testValues
#percLogevents_=Step3$percLogevents
#
#sppCumCatch=lapply(colnames(mean_),function(x) {
#                                  cum_ <- cumsum(sort(mean_[,x],decreasing=TRUE)) 
#                                  res <- cum_[seq(length(which(cum_<=thresholdCatch))+1)]
#                                  }) 
#
#sppTestValue=lapply(colnames(testValues_),function(x) which(testValues_[,x]>=thresholdTestValue))                                                             
#
#sppLogevents=lapply(rownames(percLogevents_),function(x) which(percLogevents_[x,]>=thresholdLogevents))
#
#sppAll=list()
#for (i in 1:length(sppCumCatch)){
#  sppAll[[i]]=intersect(names(unlist(sppCumCatch[[i]])),names(unlist(sppTestValue[[i]])))
#  sppAll[[i]]=intersect(sppAll[[i]],names(unlist(sppLogevents[[i]])))
#}


# Selection des espèces cibles dans les espèces présentes dans tabClusters
thresholdCatch=75
thresholdTestValue=50
thresholdLogevents=50
nbClust=Step3$nbClust
tabClusters=Step3$tabClusters


sppCumCatch=list()
sppTestValue=list()
sppLogevents=list()
sppAll=list()

for (i in 1:nbClust){
  percCatchCum=cumsum(as.numeric(tabClusters[,"% Catch",i]))
  nbSpSel=length(which(percCatchCum<thresholdCatch))+1
  sppCumCatch[[i]]=tabClusters[1:nbSpSel,"FAO",i]
  
  sppTestValue[[i]]=tabClusters[which(as.numeric(tabClusters[,"Test-value",i])>thresholdTestValue),"FAO",i]
  
  sppLogevents[[i]]=tabClusters[which(as.numeric(tabClusters[,"% Logevents",i])>thresholdLogevents),"FAO",i]
}

for (i in 1:nbClust){
  sppAll[[i]]=intersect(sppCumCatch[[i]],sppTestValue[[i]])
  sppAll[[i]]=intersect(sppAll[[i]],sppLogevents[[i]])
}
sppAll