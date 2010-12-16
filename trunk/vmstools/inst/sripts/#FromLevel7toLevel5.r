# from step 7 to step 5

path <- "C:/CLARA/VMSTools_sandbox/"
setwd(path)

#extracting the step 3

load("DK_OTB2007_3a4/PCA_70/CLARA/Step3.Rdata")

# characteristing the main species

#checking the most robust metrics

par(mfrow=c(2,2))

for (metrics in c("1st Qu.","Mean","Median","3rd Qu.")) {
mean_ <-Step3$summaryClusters[metrics,,]

#plotting the number of main species by percentage threshold
threshold <-seq(50,5,by=-5)
nb_spp <- matrix(0,ncol=Step3$nbClust,nrow=length(threshold),dimnames=list(threshold,1:Step3$nbClust))
for (i in threshold) {
spp <-lapply(colnames(mean_),function(x) which(mean_[,x]>=i))
nb_spp[as.character(i),] <- sapply(spp, length)
}
barplot(t(nb_spp),beside=T, ylab="Nb Main Species",xlab=" min % threshold in cluster",col=rainbow(Step3$nbClust),main=metrics)

}

# final choice - Median 10%
mean_ <-Step3$summaryClusters["Median",,]
spp <-lapply(colnames(mean_),function(x) which(mean_[,x]>=10))  # this gives the name of the main species in the cluster

# linking to level 5

level5 <- lapply(spp, Level5)


# Determine the metier of each logevent from dat thanks to the first species in catch of the logevent and compare with the results of classif_step3 method.
first=firstSpecies(dat1,Step2,Step3$clusters$clustering)
# Determine the metier in level 5 from the metier in level 7 based on the first species in catch (logevent by logevent)
logSpeciesL5=sapply(first$firstSp,function(x) Level7to5(x))
unique(logSpeciesL5)
barplot(table(logSpeciesL5), main="Number of logevents by group of species level 5", xlab="Groups of species level 5", ylab="Number of logevents")


# Determine the metier of each logevent thanks to the first group of species (level 5) in catch of the logevent

dat=as.matrix(dat1[,2:ncol(dat1)])
colnames(dat)=sapply(as.list(colnames(dat)), Level7to5)
nbC5=length(unique(colnames(dat)))
tab=matrix(NA,nrow=nrow(dat),ncol=nbC5)
colnames(tab)=unique(colnames(dat))
for(i in 1:nbC5){
  subdati=dat[,which(colnames(dat)==(unique(colnames(dat))[i]))]
  if(length(which(colnames(dat)==(unique(colnames(dat))[i])))==1) tab[,i]=subdati
  else
    tab[,i]=apply(subdati,1,sum)
}
# Metiers level 5 of each logevent found thanks to the first group of species in catch (level 5)
logGroupL5=unique(colnames(dat))[apply(tab,1,which.max)]

# Metiers level 5 of each logevent found thanks to the classif_step3 method.
clustersStep3L5=unlist(sapply(Step3$clusters$clustering,function(x) level5[x]))


# Compare the metiers (level 5) found thanks to the classif (Step3) vs the first group of species (level 5) in catch
#compL5ClustStep3VsGroup=table(clustersStep3L5,logGroupL5)  !! pb regroupement des clusters qui correspondent au même groupe d'espèce niveau 5.
clustStep3L5=Step3$clusters$clustering
compL5ClustStep3VsGroup=table(clustStep3L5,logGroupL5)
rownames(compL5ClustStep3VsGroup)=unlist(level5)

# Compare the metiers (level 5) found thanks to the first species in catch vs the first group of species in catch
compL5SpeciesVsGroup=table(logSpeciesL5,logGroupL5)

# Compare the metiers (level 5) found thanks to the classif (Step3) vs the first species in catch
#compL5ClustStep3VsSpecies=table(clustersStep3L5,logSpeciesL5)    !! pb regroupement des clusters qui correspondent au même groupe d'espèce niveau 5.
clustStep3L5=Step3$clusters$clustering
compL5ClustStep3VsSpecies=table(clustStep3L5,logSpeciesL5)
rownames(compL5ClustStep3VsSpecies)=unlist(level5)

