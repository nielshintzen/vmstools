################################################################################
# Calculate the cluster's within-variance (square distance between each row of #
# the cluster and the cluster's center of gravity)                             #
# we calculate the distance row by row                                         #
################################################################################

withinVar <- function(oneRowOfCluster,centerOfGravityClusti){
  comb   <- rbind(centerOfGravityClusti, oneRowOfCluster)
  sqrComb <- dist(comb)^2
return(sqrComb)}