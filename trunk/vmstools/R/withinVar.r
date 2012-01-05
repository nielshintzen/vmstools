################################################################################
# Calculate the cluster's within-variance (square distance between each row of #
# the cluster and the cluster's center of gravity)                             #
# we calculate the distance row by row                                         #
################################################################################

withinVar=function(oneRowOfCluster,centerOfGravityClusti){
  CoGWithOneRow=numeric()
  CoGWithOneRow=rbind(centerOfGravityClusti, oneRowOfCluster)
  distRowCoG=dist(CoGWithOneRow)^2
  return(distRowCoG)
}
