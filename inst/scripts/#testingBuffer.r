rlon      <- range(c(xs,outpointr[,1],outpointl[,1]),na.rm=T)
rlat      <- range(c(ys,outpointr[,2],outpointl[,2]),na.rm=T)

steplon   <- min(distance(xs[1:99],ys[1],xs[2:100],ys[1]))
steplat   <- min(distance(xs[1],ys[1:99],xs[1],ys[2:100]))

keep    <- numeric()
for(i in seq(rlon[1],rlon[2],steplon)){
  for(j in seq(rlat[1],rlat[2],steplat)){
  res <- sort.int(sqrt(abs(outer(xs,i,"-"))^2+abs(outer(ys,j,"-"))^2),index.return=T)
  points(xs[res$ix[which(res$x < sqrt(steplon^2 + steplat^2))]],ys[res$ix[which(res$x < sqrt(steplon^2 + steplat^2))]],pch=19,col="blue")